function [resSol, wellSol] = solveBlackOilWellSystem(resSol, wellSol,   ...
                                                     G, rock, S, fluid, ...
                                                     p0, dt, varargin)
%One iteration of succsubst algorithm on Black Oil pressure system.
%
% SYNOPSIS:
%   [xr, xw] = solveBlackOilWellSystem(xr, xw, G, rock, S, fluid, p0, dt)
%   [xr, xw] = solveBlackOilWellSystem(xr, xw, G, rock, S, fluid, p0, dt, ...
%                                      'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   xr, xw  - Reservoir and well solution structures either properly
%             initialized or the results from a previous call to function
%             'solveBlackOilWellSystem' (i.e., previous iteration of
%             successive substitution algorithm or previous time step),
%             and, possibly, a transport solver.
%
%   G       - Grid structure as described in grid_structure.
%
%   rock    - Rock data structure.  Must contain valid field 'rock.perm'.
%
%   S       - Linear system structure as defined by function
%             'computeMimeticIP'.
%
%   fluid   - Black Oil fluid object as defined by, e.g., function
%             'initBlackoilFluid'.
%
%   p0      - Vector, length G.cells.num, of cell pressures at previous
%             time step (not previous iteration of successive substitution
%             algorithm).
%
%   dt      - Time step size.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells    - Well structure as defined by function 'addWell'.  May be
%              empty (i.e., W = [], default value) which is interpreted as
%              a model without any wells.
%
%   bc       - Boundary condition structure as defined by function 'addBC'.
%              This structure accounts for all external boundary conditions
%              to the reservoir flow.  May be empty (i.e., bc = [], default
%              value) which is interpreted as all external no-flow
%              (homogeneous Neumann) conditions.
%
%   src      - Explicit source contributions as defined by function
%              'addSource'.  May be empty (i.e., src = [], default value)
%              which is interpreted as a reservoir model without explicit
%              sources.
%
%   LinSolve - Handle to linear system solver software to which the fully
%              assembled system of linear equations will be passed.
%              Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%              in order to solve a system Ax=b of linear equations.
%              Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%          - cellFlux     -- Outgoing flux across each local interface for
%                            all cells in the model.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% SEE ALSO:
%   initBlackoilFluid, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-07 19:22:40 +0200 (sø, 07 jun 2009) $
% $Revision: 2342 $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve', @mldivide);
   opt = merge_options(opt, varargin{:});

   % Compute gravity terms.
   grav   = gravity();
   grav   = reshape(grav(1 : size(G.nodes.coords,2)), 1, []);

   nc     = G.cells.num;
   cf     = G.cellFaces(:,1);
   cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';
   sgn    = 2*double(G.faces.neighbors(cf, 1) == cellNo) - 1;

   % cfn  = cell-face normal = *outward* face normals on each cell.
   cfn    = bsxfun(@times, G.faces.normals(cf,:), sgn);

   % ng  == n' *     g for all cellfaces.
   % nKg == n' * K * g for all cellfaces.
   [K, row, col] = permTensor(rock, size(G.nodes.coords,2));
   ng     = cfn * reshape(grav, [], 1);
   nKg    = sum(cfn(:,row) .* bsxfun(@times, K(cellNo,:), grav(col)), 2);

   %-----------------------------------------------------------------------
   % Get updated pressure and saturation dependent quantities -------------
   %
   [c, rho, mu, u] = fluid.pvt(resSol.cellPressure, resSol.z);
   s      = bsxfun(@rdivide, u, sum(u,2));
   mob    = bsxfun(@rdivide, fluid.relperm(s), mu);

   % Parameters used in pressure equation.
   g2     = grav * grav.';      % Recall: 'grav' is a *ROW* vector.
   Lt     = sum(       mob                              , 2)      ;
   ct     = sum( s  .*  c                               , 2)      ;
   omega  = sum(rho .* mob                              , 2) ./ Lt;
   zeta   = sum( c  .* mob                              , 2) ./ Lt;
   beta   = sum( c  .* mob .* bsxfun(@minus, omega, rho), 2)      ;

   %% Mimetic system with bc
   BIV    = bsxfun(@times, ng, beta(cellNo)) + ...
            zeta(cellNo) .* (resSol.cellFlux + omega(cellNo).*nKg);
   BIV    = sparse(1 : numel(cellNo), cellNo, BIV, numel(cellNo), nc);

   accum  = ct .* poreVolume(G, rock) ./ dt;
   P      = sparse(1 : nc, 1 : nc, -accum);  % == spdiags(-accum,0,nc,nc)
   g_comp = p0 .* accum;

   g_grav = - G.cells.volumes .* beta .* omega .* g2;

   %% Wells

   % Assume empty well structure by default...
   [BIW, BIVW, CW, DW, fW, hW] = deal([]);
   dFw = logical([]);
   dCw = [];

   if ~isempty(opt.wells),
      % but fill in values when there nevertheless are some...
      nperf = cellfun(@numel, { opt.wells.cells });
      wc    = vertcat(opt.wells.cells);   n = numel(wc);   i = 1 : n;
      nW    = numel(opt.wells);

      % Diagonal transmissibility matrix (inner product).
      BIW   = sparse(i, i, Lt(wc) .* vertcat(opt.wells.WI), n, n);

      % Connection matrices for all wells.
      CW    = sparse(i, wc                      , 1, n, nc);
      DW    = sparse(i, rldecode(1:nW, nperf, 2), 1, n, nW);

      % Which (and what) are the prescribed well bottom-hole pressures?
      dFw   = strcmpi('bhp', { opt.wells.type } .');
      dCw   = reshape([ opt.wells(dFw).val ], [], 1);

      % Form linsys rhs contributions, fW -> pressure, hW -> rate.
      fW    = -vertcat(opt.wells.val);  fW(~dFw) = 0;  % Remove rates
      hW    = -vertcat(opt.wells.val);  hW( dFw) = 0;  % Remove pressures

      % Expand well pressure rhs to each perforation, adjust for gravity
      % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
      dp    = norm(gravity()) * vertcat(opt.wells.dZ) .* omega(wc);
      fW    = rldecode(fW, nperf) - dp;

      BIVW  = sparse(i, wc, -zeta(wc) .* vertcat(wellSol.flux), n, nc);
      clear i n
   end

   %-----------------------------------------------------------------------
   % Build hybrid system components ---------------------------------------
   %
   n   = numel(cellNo);         i = 1 : n;
   C   = sparse(i, cellNo, 1);
   D   = sparse(i, double(G.cellFaces(:,1)), 1, n, G.faces.num);

   BI  = blkdiag(sparse(i, i, Lt(cellNo), n, n) * S.BI, BIW);     clear i n
   BIV = [       BIV ; BIVW];
   C   = [        C  ;  CW ];
   D   = blkdiag( D  ,  DW );

   [f, g, h, gp, dF, dC] = computePressureRHS(G, omega, opt.bc, opt.src);
   g = g + g_comp + g_grav + S.RHS.volume_discrepancy;

   f = [f + gp; fW];
   h = [h     ; hW];

   %-----------------------------------------------------------------------
   % Solve global hybrid system (unsymmetric) -----------------------------
   %
   dF = [dF; dFw];
   dC = [dC; dCw];

   lam     = zeros(size(dF));
   lam(dF) = dC;
   [flux, p, lam(~dF)] = schurComplement(BI, C, D(:,~dF), BIV, P, ...
                                          f, g, h(  ~dF),         ...
                                         'LinSolve', opt.LinSolve);

   %-----------------------------------------------------------------------
   % Package solution in accessible form ----------------------------------
   %
   resSol.cellPressure(:) = p;
   resSol.facePressure(:) = lam (1 : G.faces.num);
   resSol.cellFlux(:)     = flux(1 : numel(cellNo));
   resSol.faceFlux(:)     = cellFlux2faceFlux(G, flux(1 : numel(cellNo)));

   if ~isempty(opt.wells),
      i_f = size(G.cellFaces, 1);
      i_p = G.faces.num;

      for k = 1 : nW,
         wellSol(k).flux     = -flux(i_f + 1 : i_f + nperf(k));
         wellSol(k).pressure =  lam (i_p + 1 : i_p +    1    );

         i_f = i_f + nperf(k);
         i_p = i_p +    1    ;
      end
   end
end
