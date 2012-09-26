function [resSol_sub, wellSol_sub, S] = solveBOResidual(cellInx, cellInxOverlap, resSol, wellSol, ...
                                                         G, rock, S, fluid, p0, dt, varargin)
% solveBOResidual  -- Version of solveBlackOilWellSystem for solving a
%                      local problem based on cellInx using mixed repr.
%
% SYNOPSIS:
%   [resSol, wellSol] = solveBOSubSystem(cellInx, resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt)
%   [resSol, wellSol] = solveBOSubSystem(cellInx, resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt,  ...
%                                               'pn1', pv1, ...)
%
% PARAMETERS:
%   cellInx - Indices to cells defining sub domain
%  
%   resSol  - Reservoir solution structure for full fine-grid from previous time step (or
%             previous iteration of successive substitution algorithm).
%
%   wellSol - Well solution structure for all wells from previous time step (or previous
%             iteration of successive substitution algorithm).
%
%   G       - Grid structure as described in SAMSIMGRID.
%
%   rock    - Rock structure.
%
%   S       - Linear system structure as defined by function
%             assembleMimeticSystem.
%
%   W       - Well linear system structure as defined by function
%             assembleWellSystem.
%
%   PVTTAB  - PVT tables.
%
%   p0      - Vector, length G.cells.num, of cell pressures at previous
%             time step.
%
%   dt      - Time step size.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - bc  -- Boundary condtion structure as defined by function
%                        'addBC'.  This structure accounts for all external
%                        boundary contributions to the reservoir flow.
%                        Default value: bc = [] meaning all external
%                        no-flow (homogeneous Neumann) conditions.
%
%               - src -- Explicit source contributions as defined by
%                        function 'addSource'.
%                        Default value: src = [] meaning no explicit
%                        sources exist in the model.
%
%               - LinSolve -- Solver for the resulting linear systems.
%                             Handle to a function supporting the syntax
%                                 x = LinSolve(A,b)
%                             for solving the system Ax=b of linear
%                             equations.
%                             Default value: LinSolve = @mldivide
%                             (backslash).
%
% RETURNS:
%   resSol  - Updated sparse reservoir solution structure.
%
%   wellSol - Updated 'sparse' well solution structure.


%--------------------------------------------------------------------------
% Get updated pressure and sat-dependent quantities -----------------------
%

opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve', @mldivide);
opt = merge_options(opt, varargin{:});
   



%-----------------------------------------------------------------------
% Get updated pressure and saturation dependent quantities -------------
[c, rho, mu, u] = fluid.pvt(resSol.cellPressure, resSol.z);
s      = bsxfun(@rdivide, u, sum(u,2));
mob    = fluid.relperm(s) ./ mu;   
Lti    = 1./sum(mob, 2);   


% System matrices
nc     = G.cells.num;
cf     = G.cellFaces(:,1);
ncf    = numel(cf);
cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';

B  = spdiags(Lti(cellNo), ncf, ncf) * S.B;
C  = sparse((1:ncf)', cellNo, 1, ncf, nc);
D  = sparse((1:ncf)', double(cf), 1, ncf, G.faces.num);
accum  = sum(c, 2) .* poreVolume(G, rock) ./ dt;
P      = spdiags(-accum,0,nc,nc)

% System RHS
grav   = gravity(); grav = grav(:);
hg     = ( G.faces.centroids(cf, :) - G.cells.centroids(cellNo, :) ) * grav;
rhoTot = sum(rho.*mob ,2).*Lti;
f      = hg .* rhoTot(cellNo);

[f, g, h, gp, dF, dC] = computePressureRHS(G, rhoTot, opt.bc, opt.src);
a1   =  sum( c .* mob, 2) .* Lti;
a2   =  sum( c  .* mob .* bsxfun(@minus, 2*rhoTot, rho), 2) .* Lti;  % sjekk fortegn???
a3   =  rhoTot .* sum( c .* mob .* bsxfun(@minus, rhoTot, rho));     % sjekk fortegn???





g = g + g_comp + g_grav + S.RHS.volume_discrepancy;

f = [f + gp; fW];
h = [h     ; hW];


[f, g, h] = linHybRHS()
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

   %--------------------------------------------------------------------------
   % residual RHS --------------------------------
   %
   v   = resSol.cellFlux;
   p   = resSol.cellPressure;
   sD  = sum(D)';
   Dp  = spdiags(1./sD, 0, numel(sD), numel(sD));

   pi = Dp*D'*(- BI\v + C*p);
   %pi  = D \ (- BI\v + C*p);
   %fr  = - BI\v + C*p - D*pi;
   gr  = g - C'*v + BIV'*(BI\v) + P*p; %(p-p0);
   hr  = h - D'*v;

   %--------------------------------------------------------------------------
   % Produce local indices --------------------------------
   %
   if ~islogical(cellInx)
      indC = logical( sparse(cellInx, 1, 1, G.cells.num, 1) );
   else
      indC = cellInx;
   end
   if ~islogical(cellInxOverlap)
      indCO = logical( sparse(cellInxOverlap, 1, 1, G.cells.num, 1) );
   else
      indCO = cellInxOverlap;
   end
   indB     = logical( C * indC );
   indD     = logical( D' * indB );

   indBO    = logical( C * indCO );
   indDO    = logical( D' * indBO );

   %--------------------------------------------------------------------------
   % Solve local hybrid system (unsymmetric) --------------------------------
   %
   %{
   if ~isempty(W),
      dFw = strcmp({ W.type } .', 'bhp');
      dCw = [ W(dFw).val ].';
   else
      dFw = logical([]);
      dCw = [];
   end
   %}
   dF = [dF; dFw];
   dC = [dC; dCw];

   % sub matrices:
   %indDD  = and(indD, dF);
   %----------------------------------
   indDnD = and(indDO, ~dF);
   %indDnD = indDO;

   fr = f - BI\v + C*p - D(:, ~dF)*pi(~dF);
   %fr = - BI\v + C*p - D*pi;
   %----------------------------------
   BIs  = BI(indBO, indBO);
   Cs   = C(indBO, indCO);
   %dFs  = dF(indD);
   Ds   = D(indBO, indDnD );
   BIVs = BIV(indBO, indCO);
   Ps   = P(indCO, indCO);
   fl   = fr.*indB; fs   = fl(indBO);
   gl   = gr.*indC; gs   = gl(indCO);
   %----------------------------------
   hs   = hr( indDnD );
   %hs   = zeros(size(hs));
   %----------------------------------

   flux    = sparse(size(S.B,2), 1);
   p       = sparse(G.cells.num, 1);
   lam     = sparse(size(D,2)  , 1);
   lam(dF) = dC;
   lam     = lam.*indD;

   [flux(indBO), p(indCO), lam(indDnD)] = schurComplement(BIs, Cs, Ds, BIVs, Ps, ...
                                                           fs, gs, hs,           ...
                                                          'LinSolve', opt.LinSolve);

   %--------------------------------------------------------------------------
   % Package solution in accessible form -------------------------------------
   %

   resSol_sub.cellPressure = p;
   resSol_sub.facePressure = lam (1 : G.faces.num);
   resSol_sub.cellFlux     = flux(1 : size(S.B,1));
   resSol_sub.faceFlux     = cellFlux2faceFlux(G, flux(1 : size(S.B,1)));

   if ~isempty(opt.wells),
      i_f = size(G.cellFaces, 1);
      i_p = G.faces.num;

      for k = 1 : nW,
         wellSol_sub(k).flux     = -flux(i_f + 1 : i_f + nperf(k)); %#ok
         wellSol_sub(k).pressure =  lam (i_p + 1 : i_p +    1    ); %#ok

         i_f = i_f + nperf(k);
         i_p = i_p +    1    ;
      end
   else
      wellSol_sub = struct([]);
   end


   
