function [xr, xw] = solveLinBO(xr, xw, g, rock, s, fluid, p0, dt, varargin)
%Solve incompressible flow problem (fluxes/pressures).
%
% SYNOPSIS:
%   [xr, xw] = solveIncompFlow(xr, xw, G, S, fluid)
%   [xr, xw] = solveIncompFlow(xr, xw, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (wells, sources, and boundary conditions).
%
% REQUIRED PARAMETERS:
%   xr, xw - Reservoir and well solution structures either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'solveIncompFlow' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, S   - Grid and (mimetic) linear system data structures as defined by
%            function 'computeMimeticIP'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   Solver       - Which solver mode function 'solveIncompFlow' should
%                  employ in assembling and solving the block system of
%                  linear equations.
%                  String.  Default value: Solver = 'hybrid'.
%
%                  Supported values are:
%                    - 'hybrid' --
%                         Assemble and solve hybrid system for interface
%                         pressures.  System is eventually solved by Schur
%                         complement reduction and back substitution.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','hybrid') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
%
%                    - 'mixed' --
%                         Assemble and solve a hybrid system for interface
%                         pressures, cell pressures and interface fluxes.
%                         System is eventually reduced to a mixed system as
%                         per function 'mixedSymm'.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','mixed') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
%
%                    - 'tpfa' --
%                         Assemble and solve a cell-centred system for cell
%                         pressures.  Interface fluxes recovered through
%                         back substitution.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','mixed') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
%
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput - Whether or not to return the final system matrix 'A' to
%                  the caller of function 'solveIncompFlow'.
%                  Logical.  Default value: MatrixOutput = FALSE.
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
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
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
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'solveIncompFlow:DrivingForce:Missing'
%
%
% SEE ALSO:
%   computeMimeticIP, addBC, addSource, addWell, initSimpleFluid
%   initResSol, initWellSol, solveIncompFlowMS.

%{
#COPYRIGHT#
%}

% $Id: solveIncompFlow.m 2342 2009-06-07 17:22:40Z bska $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'Solver',       'mixed',         ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false);
   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   no_grav = ~(norm(g_vec(1 : size(g.nodes.coords,2))) > 0);
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), no_grav]),
      warning(id('DrivingForce:Missing'),                      ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   else
      cellNo  = rldecode(1:g.cells.num, double(g.cells.numFaces), 2) .';
      s.C     = sparse(1:numel(cellNo), cellNo, 1);
      s.D     = sparse(1:numel(cellNo), double(g.cellFaces(:,1)), 1, ...
                       numel(cellNo), g.faces.num);
      s.sizeB = repmat(size(g.cellFaces, 1), [1,2]);
      s.sizeC = size(s.C);
      s.sizeD = size(s.D);

      [A, b, dF, dC] = build_system(xr, xw, g, s, rock, fluid, p0, dt, ...
                                    opt.wells, opt.bc, opt.src, opt);

      solver   = pick_solver(g, s, dF, dC, opt);
      x        = solver(A, b);

      [xr, xw] = pack_solution(xr, xw, g, s, x{1:3}, opt);
      if opt.MatrixOutput, xr.A = x{4}; end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function s = id(s)
   s = ['solveIncompFlow:', s];
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = build_system(xr, xw, g, s, rock, fluid, p0, dt, ...
                                       w, bc, src, opt)
   % Build system components B, C, D, f, g, and h for final hybrid system
   %
   %    [  B   C   D  ] [  v ]     [ f ]
   %    [  C'  P   0  ] [ -p ]  =  [ g ] .
   %    [  D'  0   0  ] [ cp ]     [ h ]
   %
   % The components are represented as one-dimensional cell arrays as
   %
   %    A = { B, C, D ,P }
   %    b = { f, g, h }
   %
   % both of which are assembled separately from reservoir contributions
   % (syscomp_res) and well contributions (syscomp_well).  We assemble the
   % global components by concatenating the well components onto the
   % reservoir components.  This assembly is a result of treating wells as
   % faces/contacts.
   %
   % Specifically, we assemble the components as
   %
   %    A{1} = B = [ Br, 0; 0, Bw ]  % BLKDIAG
   %    A{2} = C = [ Cr   ;    Cw ]  % VERTCAT
   %    A{3} = D = [ Dr, 0; 0, Dw ]  % BLKDIAG
   %    A{4} = P = Pr
   %
   %    b{1} = f = [ fr   ;    fw ]  % VERTCAT
   %    b{2} = g =   gr + gw         % (sparse contribution from gw)
   %    b{3} = h = [ hr   ;    hw ]  % VERTCAT
   %
   [Ar, br, dFr, dCr] = syscomp_res  (xr, g, s, rock, fluid, p0, dt, ...
                                      bc, src, opt);
   [Aw, bw, dFw, dCw] = syscomp_wells(xr, xw, g, s, w, rock, fluid, opt);

   A    = cell([1, 4]);           b    = cell([1, 3]);
   A{1} = blkdiag(Ar{1}, Aw{1});  b{1} = vertcat(br{1}, bw{1});
   A{2} = vertcat(Ar{2}, Aw{2});  b{2} = br{2} + full(bw{2})  ;
   A{3} = blkdiag(Ar{3}, Aw{3});  b{3} = vertcat(br{3}, bw{2});
   A{4} = Ar{4};

   % Finally, assemble a list of prescribed (contact) pressures, i.e., both
   % known face pressure and known bottom hole pressures.  These will be
   % eliminated from the linear system prior to calling one of the system
   % solvers, and entered directly into the solution component 'cp'.
   %
   dF   = [dFr; dFw];
   dC   = [dCr; dCw];
end

%--------------------------------------------------------------------------

function solver = pick_solver(g, s, dF, dC, opt)
   regul = ~any(dF);  % Set zero level if no prescribed pressure values.

   switch lower(opt.Solver),
      case 'hybrid',
          error(id('SolverMode:Inconsistent'), ...
                 ['solveLinBO only works with system type', ...
                  ' ''mixed'' for now'], s.type);
         
         %if ~any(strcmp(s.type, {'hybrid', 'comp_hybrid'})),
         %   error(id('SolverMode:Inconsistent'), ...
         %        ['Solver mode ''hybrid'' is incompatible with ', ...
         %         'linear system type ''%s''.'], s.type);
         %end

         %solver = @solve_hybrid;
      case 'mixed',
         if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                 ['Solver mode ''mixed'' is incompatible with ', ...
                  'linear system type ''%s''.'], s.type);
         end

         solver = @(A,b) solve_mixed(A, b, @solveMixedLinSys);
      case 'tpfa',
          error(id('SolverMode:Inconsistent'), ...
                 ['solveLinBO only works with system type', ...
                  ' ''mixed'' for now'], s.type);
         %if ~any(strcmp(s.type, {'mixed', 'tpfa', 'comp_hybrid'})),
         %   error(id('SolverMode:Inconsistent'), ...
         %        ['Solver mode ''tpfa'' is incompatible with ', ...
         %         'linear system type ''%s''.'], s.type);
         %end

         %solver = @(A,b) solve_mixed(A, b, @tpfSymm);
      otherwise,
         error(id('SolverMode:NotSupported'), ...
               'Solver mode ''%s'' is not supported.', opt.solver);
   end

   % Specify prescribed pressures on 'Dirichlet' contacts.
   lam     = zeros(size(dF));
   lam(dF) = dC;

   function x = solve_hybrid(A, b)
      A{3} = A{3}(:,~dF);
      b{3} = b{3}(  ~dF);
      [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                     'Regularize', regul, ...
                                     'LinSolve', opt.LinSolve);
      lam(~dF) = x{3};
      x{3}     = lam;
   end

   function x = solve_mixed(A, b, solver)
      nF = neumann_faces(g, dF, opt);
      cellNo = rldecode((1 : g.cells.num) .', double(g.cells.numFaces));
      sgn = 2*double(g.faces.neighbors(g.cellFaces(:,1), 1) == cellNo) - 1;

      Do = oriented_mapping(s, sgn, opt);

      A{3} = A{3}(:,nF);
      b{3} = b{3}(  nF);
      [x{1:4}] = solver(A{:}, b{:}, Do, 'Regularize', regul, ...
                                  'LinSolve', opt.LinSolve);

      x{1}    = [faceFlux2cellFlux(g, x{1}(1 : g.faces.num)); ...
                 x{1}(g.faces.num + 1 : end)];

      lam(nF) = x{3};
      x{3}    = lam;
   end
end

%--------------------------------------------------------------------------

function [xr, xw] = pack_solution(xr, xw, g, s, flux, pres, lam, opt)
   xr.cellPressure(:) = pres;
   xr.facePressure(:) = lam (1 : s.sizeD(2));
   xr.cellFlux(:)     = flux(1 : s.sizeB(1));
   xr.faceFlux(:)     = cellFlux2faceFlux(g, flux(1 : s.sizeB(1)));

   if ~isempty(opt.wells),
      nw  = numel(opt.wells);
      i_f = s.sizeB(1);
      i_p = s.sizeD(2);

      for k = 1 : nw,
         nperf = numel(opt.wells(k).cells);

         xw(k).flux     = -flux(i_f + 1 : i_f + nperf);
         xw(k).pressure =  lam (i_p + 1 : i_p +   1  );

         i_f = i_f + nperf;
         i_p = i_p +   1  ;
      end
   end
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_res(xr, G, S, rock, fluid, p0, dt, ...
                                      bc, src, opt)
   A = cell([1, 4]);

   % Get updated pressure and saturation dependent quantities -------------
   [c, rho, mu, u] = fluid.pvt(xr.cellPressure, xr.z);
   s      = bsxfun(@rdivide, u, sum(u,2));
   ct     = sum( s  .*  c , 2);
   mob    = fluid.relperm(s) ./ mu;
   mobt   = sum(mob, 2);
   rhoTot = sum(rho.*mob ,2)./mobt;
   
   % System matrices ---------------------------------------------------------
   nc     = G.cells.num;
   cf     = G.cellFaces(:, 1);
   ncf    = numel(cf);
   
   Mobt = spdiags(S.C * mobt, 0, S.sizeB(1), S.sizeB(2));
   
   if strcmpi(opt.Solver, 'hybrid'),
       if ~isfield(S, 'BI'),
           error(id('SolverMode:Inconsistent'), ...
               ['Solver mode ''hybrid'' is incompatible with ', ...
               'linear system type ''%s''.'], s.type);
       end
       A{1} = Mobt * S.BI;
   else
       if ~isfield(S, 'B'),
           error(id('SolverMode:Inconsistent'), ...
               ['Solver mode ''%s'' is incompatible with ', ...
               'linear system type ''%s''.'], opt.Solver, s.type);
       end
       A{1} = Mobt \ S.B;
   end
   
   A{2} = S.C;
   A{3} = S.D;
   accum   = ct .* poreVolume(G, rock) ./ dt;
   A{4}    = spdiags(-accum, 0, G.cells.num, G.cells.num);

   % System RHS  ----------------------------------------------------------
   grav = gravity(); grav = grav(:);
   [f, g, h, fg, dF, dC] = computePressureRHS(G, rhoTot, bc, src);
   f = f + fg;
   %[b{1:3}, grav, dF, dC] = computePressureRHS(g, omega, bc, src);
   %b{1} = b{1} + grav;  % Add (mimetic) gravity contributions.
   % Accumulation
   g = g + accum .* p0;
   
   % a1 - part
   a1   =  sum( c .* mob, 2) ./ mobt;
   %vv   = sparse((1:ncf)',cellNo, xr.cellFlux, ncf, nc);
   %g    = g + a1 .* ( vv' * ( B*xr.cellFlux ) );
   g    = g + a1 .* ( S.C' * ( xr.cellFlux .* (A{1}*xr.cellFlux) ) );
   
   if any(grav)
       % a2 - part
       hg = ( G.faces.centroids(cf, :) - G.cells.centroids(cellNo, :) ) * grav;
       a2 =  sum( c  .* mob .* bsxfun(@minus, 2*rhoTot, rho), 2) ./ mobt;  % sjekk fortegn???
       g  = g - a2 .* (vv'*hg);
       % a3 - part
       a3   =  rhoTot .* sum( c .* mob .* bsxfun(@minus, rhoTot, rho));     % sjekk fortegn???
       [K, row, col] = permTensor(rock, size(G.nodes.coords,2));
       gKg  = bsxfun(@times, K, grav(col)) .* grav(row)';
       g = g + a3 .* poreVolume(G, rock) .* gKg;
   end
   b{1} = f; b{2} = g; b{3} = h;
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_wells(xr, xw, G, S, W, rock, fluid, opt)
   % Assume empty well structure by default...
   A  = cell([1, 3]);
   b  = cell([1, 3]);
   b{2} = sparse(G.cells.num, 1);
   dF = logical([]);
   dC = [];

   if ~isempty(W),
       nperf = cellfun(@numel, { W.cells });
       wc    = vertcat(W.cells);   np = numel(wc);
       nw    = numel(W);
       rates = vertcat( xw.flux );
       
       [c, rho, mu, u] = fluid.pvt(xr.cellPressure(wc), xr.z(wc));
       s      = bsxfun(@rdivide, u, sum(u,2));
       mob    = fluid.relperm(s(wc)) ./ mu;
       mobt   = sum(mob, 2);
       rhoTot = sum(rho.*mob ,2)./mobt;
       a1     = sum( c .* mob, 2) ./mobt;
       
       % Diagonal transmissibility matrix (inner product).
       A{1}   = sparse(1:np, 1:np, (1./mobt) .* 1./( vertcat(W.WI) ), np, np);
       
       % Connection matrices for all wells.
       A{2}   = sparse(1:np, wc , 1, np, G.cells.num);
       A{3}   = sparse(1:np, rldecode(1:nw, nperf, 2), 1, np, nW);
       
       % Which (and what) are the prescribed well bottom-hole pressures?
       dF   = strcmpi('bhp', { W.type } .');
       
       % Form linsys rhs contributions, fW -> pressure, hW -> rate.
       f = -vertcat(W.val);  f(~dF) = 0;  % Remove rates
       h = -vertcat(W.val);  h( dF) = 0;  % Remove pressures
       g = sparse(G.cells.num, 1);
       g(wc) = a1 .* rates .* (A{1} * rates);
       
       % Expand well pressure rhs to each perforation, adjust for gravity
       % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
       dp  = computeDp(W, rho, rhoTot, wellSol);
       f   = rldecode(f, nperf) - dp;
       b{1} = f; b{2} = g; b{3} = h;
       %{
      % but fill in values when there nevertheless are some...
      nperf = cellfun(@numel, { W.cells });
      wc    = vertcat(W.cells);   n = numel(wc);   i = 1 : n;
      nW    = numel(W);

      % Diagonal transmissibility matrix (inner product).
      v = mob(wc) .* vertcat(W.WI);
      if ~strcmpi(opt.Solver, 'hybrid'), v = 1 ./ v; end
      A{1} = sparse(i, i, v, n, n); % == spdiags(v,0,n,n), but more direct.

      % Connection matrices {2} -> C and {3} -> D for all wells.
      A{2} = sparse(i, wc                      , 1, n, G.cells.num);
      A{3} = sparse(i, rldecode(1:nW, nperf, 2), 1, n, nW         );

      % Which (and what) are the prescribed well bottom-hole pressures?
      dF   = strcmpi('bhp', { W.type } .');
      dC   = reshape([ W(dF).val ], [], 1);

      % Form linsys rhs contributions, {1} -> pressure, {2} -> rate.
      b{1} = -vertcat(W.val);  b{1}(~dF) = 0;  % Remove rates
      b{2} = -vertcat(W.val);  b{2}( dF) = 0;  % Remove pressures

      % Expand well pressure rhs to each perforation, adjust for gravity
      % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
      dp   = norm(gravity()) * vertcat(W.dZ) .* omega(wc);
      b{1} = rldecode(b{1}, nperf) - dp;
       %}
   end
end

%--------------------------------------------------------------------------

function nF = neumann_faces(g, dF, opt)
% Determine the 'Neumann faces' (nF) (i.e., the faces for which
% face/contact pressures will be determined in the 'mixed' and 'TPFA'
% cases) of the model (g) according to the following rules:
%
%   - No internal faces are counted amongst the 'Neumann faces'.
%   - An external face is a 'Neumann face' unless there is a prescribed
%     pressure value (i.e., a Dirichlet condition) associated to the face.
%   - Wells controlled by rate-constraints are treated as Neumann contacts.
%
   nF                                 = false([g.faces.num, 1]); % No int.
   nF(any(g.faces.neighbors == 0, 2)) = true;  % All external faces ...
   nF(dF(1 : g.faces.num))            = false; % ... except Dirichlet cond.

   if ~isempty(opt.wells),
      % Additionally include all rate-constrained wells.
      nF = [nF; strcmpi('rate', { opt.wells.type } .')];
   end
end

%--------------------------------------------------------------------------

function Do = oriented_mapping(s, orient, opt)
% 'Do' maps face fluxes to half-face fluxes.  This matrix is used to form
% the reduced, mixed system of linear equtions in linear system solver
% functions 'mixedSymm' and 'tpfSymm'.
%
   if ~isempty(opt.wells),
      n  = sum(cellfun(@numel, { opt.wells.cells }));
      dw = speye(n);
      ow = ones([n, 1]);
   else
      dw = [];
      ow = [];
   end
   nf = numel(orient) + numel(ow);
   Do = spdiags([orient; ow], 0, nf, nf) * blkdiag(s.D, dw);
end

function dp = computeDp(W, rho, rhoTot, wellSol)
% Find pressure drop distribution along well. Assume all phases flow freely
% at infinite speed such that segment volume fraction of phase a is equal to
% inflow fraction of phase a. Pressure drop in a segment is then rho*g*Dz,
% where rho is the mixture density. Further for now make simplifying assumption
% that densities in segment equal densities in well-cell.
nperf = cellfun(@numel, { W.cells });
dp    = zeros(prod(nperf), 1);
grav  = gravity();
if norm(grav)>0
    inx = 0;
    for k = 1 : numel(W)
        wc = W(k).cells;
        numSeg   = numel(wc);
        inxs = inx + (1:nperf(k))';
        inx  = inx + nperf(k);
        if numSeg > 1
            resFlux  = wellSol(k).flux;
            totFlux  = sum( resFlux );
            if any(sign(resFlux) == -sign(totFlux)  )
                warning('Both in- and out-flux in same well.')
            end
            wellFlux = [totFlux; totFlux - cumsum( resFlux(1:end) )];
            
            rhoRes  = rhoTot(inxs);
            rhoInj  = rho(inxs(1),:) * W(k).Comp_i';
            
            res2segFlux    = resFlux .* (resFlux < 0);
            seg2segPosFlux = wellFlux .* (wellFlux > 0); %heel to toe (i -> i+1)
            seg2segNegFlux = wellFlux .* (wellFlux < 0); %toe to heel (i -> i-1)
            
            rhoSeg = rhoRes;
            
            % iterate until fixed
            rhoSegPrev = rhoSeg*0;
            while rhoSegPrev ~= rhoSeg
                rhoSegPrev = rhoSeg;
                rhoSeg = res2segFlux .* rhoRes + ...
                    seg2segPosFlux(1:numSeg) .* [rhoInj; rhoSeg(1:numSeg-1)] +...
                    seg2segNegFlux(2:numSeg+1) .* [rhoSeg(2:numSeg); 0];
            end
            perfDZ = W(k).dZ(2:numSeg) - W(k).dZ(1:numSeg-1);
            perfDp = norm(gravity()) * ...
                [0; .5 * perfDZ * (rhoSeg(1:numSeg-1) + rhoSeg(2:numSeg))];
            dp(inxs) = cumsum( perfDp );
        end
    end
end
end