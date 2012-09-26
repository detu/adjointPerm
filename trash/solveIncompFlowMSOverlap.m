function [xr, xw] = solveIncompFlowMSOverlap(xr, xw, G, CG, p, S, CS, ...
                                             fluid, varargin)
% solveIncompFlowMSOverlap -- Solve coarse (multiscale) pressure system.
%
% SYNOPSIS:
%  [xr, xw] = solveIncompFlowMSOverlap(xr, xw, G, CG, p, S, CS, W, fluid)
%  [xr, xw] = solveIncompFlowMSOverlap(xr, xw, G, CG, p, S, CS, W, ... 
%                                      fluid, 'pn1', pv1, ...)
%
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
%            function 'twophaseUpwReorder'.
%
%   G, CG, p - Grid, coarse grid, and cell-to-block partition vector,
%              respectively.
%
%   S, CS  - Linear system structure on fine grid (S) and coarse grid
%            (CS) as defined by functions 'assembleMimeticSystem' and
%            'generateCoarseSystem', respectively.
%
%   fluid  - Fluid object as defined by function initSimpleFluid.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - well linear system structure (W) as defined by functions 
%            'addWell', 'assembleWellSystem', and
%            'assembleCoarseWellSystem'.  An empty well linear system
%            structure (i.e., W=struct([])) is supported.  The resulting
%            model is interpreted as not containing any well driving
%            forces.
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
%                         'assembleMimeticSystem'.
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
%                         'assembleMimeticSystem'.
%
%                    - 'tpfa' --
%                         Assemble and solve a cell-centred system for cell
%                         pressures.  Interface fluxes recovered through
%                         back substitution.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','mixed') or option
%                         pair ('Type','comp_hybrid') to function
%                         'assembleMimeticSystem'.
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
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'. Except for
%                            mixed..
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
%           'solveIncompFlowMSOverlap:DrivingForce:Missing'
%
%
% SEE ALSO:
%   assembleMimeticSystem, addBC, addSource, addWell, initSimpleFluid
%   initResSol, initWellSol, solveIncompFlow.

% $Id: solveIncompFlowMS.m 1293 2009-02-12 14:40:08Z ilig $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'Solver',       'hybrid',         ...
                'LinSolve',     @mldivide, ...
                'MatrixOutput', false);
   opt = merge_options(opt, varargin{:});

    if all([isempty(opt.bc)   , ...
            isempty(opt.src)  , ...
            isempty(opt.wells), gravity() == 0]),
      warning(id('DrivingForce:Missing'),                       ...
              ['No external driving forces present in model--', ...
               'state remains unchanged.\n']);
    else

       [A, b, dF, dC, Phi, Psi] = build_system(xr, G, CG, p, S, CS, ...
                                               fluid, opt);
       solver   = pick_solver(CG, S, CS, dF, dC, opt);
       x        = solver(A, b);
       [xr, xw] = pack_solution(xr, xw, G, CG, p, S, CS, x{1:3}, ...
                                Phi, Psi, fluid, opt);

       if opt.MatrixOutput, xr.A = x{4}; end
    end
end
    
    
function [A, b, dF, dC, Phi, Psi] = build_system(xr, G, CG, p, S, CS, ...
                                                 fluid, opt)
   W = opt.wells;
   [nsub, sub] = get_subfaces(G, CG, CS);
   [I, J]      = coarsening_ops(G, CG, CS, sub, nsub);

   %-----------------------------------------------------------------------
   % Assemble well sub matrices -------------------------------------------
   %
   Lt = fluid.Lt(xr);
   [BW, Psi_w, Phi_w, CW, DW, fW, hW] = ...
      unpackWellSystemComponentsMS(W, G, p, Lt);
 %-----------------------------------------------------------------------
   % Eliminate known contact and well pressures from linsys ---------------
   %
   if ~isempty(W),
      dFw = reshape(strcmp({ W.type } .', 'bhp'), [], 1);
      dCw = reshape([ W(dFw).val ]              , [], 1);
   else
      dFw = logical([]);
      dCw = [];
   end   
   
   omega = fluid.omega(xr);

   %-----------------------------------------------------------------------
   %% Build full (block) system of linear equations -----------------------
   %
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;
   
   Phi     = res_basisP(G, CG, CS, p, Lt);
  
   %-----------------------------------------------------------------------
   %% Assemble (inverse) coarse mass matrix -------------------------------
   %-----------------------------------------------------------------------
   % Original mass matrix Bc = psi.' * (Bf/mob) * psi ---------------------
   %
   
   
   % Assemble linear system - Add well contributions to linear system.
   A    = cell([1, 3]);                  b    = cell([1, 3]);
   C = vertcat(CS.C(actB,:),    CW{:});
   
   
   
   if strcmp(opt.Solver, 'hybrid')
           flx_bas = [CS.basis(:, actB), -Psi_w];
      Lti = spdiags(S.C * (1 ./ Lt), 0, S.sizeB(1), S.sizeB(2));
      Psi     = S.BI * flx_bas;
      
      B   = (Lti * flx_bas).' * Psi;

      wellRange = numel(actB) + 1 : size(B,1);
      B(wellRange, wellRange) = B(wellRange, wellRange) + blkdiag(BW{:});

      %--------------------------------------------------------------------
      % Construct BI = INV(B) which is needed when solving hybrid system --
      %
      BI = invert(B, C);
      A{1} = BI;
      
      [f, g, h, dFr, dCr] = sys_rhs(G, S, omega, Psi(:, 1 : numel(actB)), ...
                                    I, J, opt.src, opt.bc);
   else


      % If basis is on (fine) faces Do extedens to (fine) cellFaces!
      %orientation = G.cellFaces(:,2);
      cellNo      = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
      orientation = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo) - 1;

      numcf       = numel(orientation);
      Do          = spdiags(orientation, 0, numcf, numcf) * S.D;

      % both reservoir basis and wellbasis on faces
      if size(Psi_w, 1) == G.faces.num
         Psi = Do* [CS.basis(:, actF), -Psi_w];
      else % well basis is on cellFaces
         Psi = [Do* CS.basis(:, actF), -Psi_w];
      end

      Lti = spdiags(S.C * (1 ./ Lt), 0, S.sizeB(1), S.sizeB(2));

      B   = Psi' * ( Lti*S.B ) * Psi;

      wellRange = numel(actF) + 1 : size(B,1);
      B(wellRange, wellRange) = B(wellRange, wellRange) + blkdiag(BW{:});

      A{1} = B;
      [f, g, h, dFr, dCr] = sys_rhs(G, S, omega, Psi(:, 1 : numel(actF)), ...
                                    I, J, opt.src, opt.bc);

      f = f(CG.cellFaces(actB,2));
   end

   
   
   % Build remainding linear system                              
                                             b{1} = vertcat(f, fW{:});
   A{2} = C;                                 b{2} = g                ;
   A{3} = blkdiag(CS.D(actB,actF), DW{:});   b{3} = vertcat(h, hW{:});
    
   dF = [dFr; dFw];
   dC = [dCr; dCw];   
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [nsub, sub] = get_subfaces(G, CG, CS)
% get_subfaces -- Extract fine-scale faces constituting all active coarse
%                 faces.
%
% SYNOPSIS:
%   [nsub, sub] = get_subfaces(G, CG, CS)
%
% PARAMETERS:
%   G, CG - Fine and coarse grid discretisation of the reservoir model,
%           respectively.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
% RETURNS:
%   nsub, sub - The return values from function 'subFaces' compacted to
%               only account for active coarse faces (i.e. coarse faces
%               across which there is a non-zero flux basis function and,
%               therefore a possibility of non-zero flux).

   [nsub, sub] = subFaces(G, CG);
   sub_ix      = cumsum([0; nsub]);
   sub         = sub(mcolon(sub_ix(CS.activeFaces) + 1, ...
                            sub_ix(CS.activeFaces + 1)));
   nsub        = nsub(CS.activeFaces);
end

%--------------------------------------------------------------------------

function [I, J] = coarsening_ops(G, CG, CS, sub, nsub)
% coarsening_ops -- Compute coarsening operators from fine to coarse model.
%
% SYNOPSIS:
%   [I, J] = coarsening_ops(G, CG, CS, sub, nsub)
%
% PARAMETERS:
%   G, CG     - Fine and coarse grid discretisation of the reservoir model,
%               respectively.
%
%   CS        - Coarse linear system structure as defined by function
%               'generateCoarseSystem'.
%
%   sub, nsub - Packed representation of fine-scale faces constituting
%               (active) coarse faces.  Assumed to follow the conventions
%               of the return values from function 'subFaces'.
%
% RETURNS:
%  I - Coarse-to-fine cell scatter operator (such that I.' is the
%      fine-to-coarse block gather (sum) operator).
%
%  J - Coarse-to-fine face scatter operator for active coarse faces such
%      that J.' is the fine-to-coarse face gather (sum) operator.

   I = CG.cells.subCells;
   J = sparse(sub, rldecode(CS.activeFaces(:), nsub), ...
              1, G.faces.num, CG.faces.num);
   J = J(:, CS.activeFaces);
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = sys_rhs(G, S, omega, Psi, I, J, src, bc)
% sys_rhs -- Evaluate coarse system right hand side contributions.
%
% SYNOPSIS:
%   [f, G, h, dF, dC] = sys_rhs(G, S, omega, Psi, I, J, src, bc)
%
% PARAMETERS:
%   G, S  - Grid and (mimetic) linear system data structures as defined by
%           function 'assembleMimeticSystem'.
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
%   Psi   - Flux basis functions from reservoir blocks.
%
%   I, J  - Coarsening operators as defined by function 'coarsening_ops'.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Coarse-scale direct (block) right hand side contributions as
%             expected by the linear system solvers such as
%             'schurComplement'.
%
%   dF, dC  - Coarse-scale packed Dirichlet/pressure condition structure.
%             The faces in 'dF' and the values in 'dC'.  May be used to
%             eliminate known face pressures from the linear system before
%             calling a system solver (e.g., 'schurComplement').
%
% SEE ALSO:
%   computePressureRHS.

   % Evaluate fine-scale system right hand side contributions.
   %
   [ff, gg, hh, dF, dC] = computePressureRHS(G, omega, bc, src);

   % Accumulate fine-scale boundary conditions to coarse scale.
   %
   ff = Psi.' * ff;
   gg =  I .' * gg;
   hh  = J .' * hh;

   % Account for coarse scale Dirichlet boundary conditions.
   %
   dC2 = zeros(size(dF));  dC2(dF) = dC;
   dF  = J.' * dF;
   dC  = J.' * dC2 ./ dF;
   dF  = logical(dF);
   dC  = dC(dF);
end

%--------------------------------------------------------------------------

function BI = invert(B, C)
   % Pre-allocate storage for sparse input arrays.
   %
   n = sum(sum(C > 0) .^ 2);
   i = zeros([n, 1]);
   j = zeros([n, 1]);
   s = zeros([n, 1]);

   % Compute inverse mass matrix values by inverting local matrices
   % (one dense matrix for each coarse block).
   %
   ix = 0;
   for k = 1 : size(C,2),
      r   = reshape(find(C(:,k)), [], 1);
      nF  = numel(r);
      nF2 = nF^2;

      i(ix + 1 : ix + nF2) = repmat(r  , [nF, 1]);
      j(ix + 1 : ix + nF2) = repmat(r.', [nF, 1]);
      s(ix + 1 : ix + nF2) = inv( B(r,r) );

      ix = ix + nF2;
   end

   % Assemble final inverse mass matrix.
   %
   BI = sparse(i, j, s, size(B,1), size(B,2));
end

%--------------------------------------------------------------------------

function Phi = res_basisP(G, CG, CS, p, mob)
% res_basisP -- Assemble pressure basis functions for reservoir blocks.
%
% SYNOPSIS:
%   Phi = res_basisP(G, CG, CS, p, mob)
%
% PARAMETERS:
%   G, CG - Fine-scale and coarse-scale discretisations (grids) of
%           reservoir model.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
%   p     - Cell-to-coarse block partition mapping vector.
%
%   mob   - Vector of (current) total mobility values.  One scalar value
%           for each cell in the underlying fine-scale model.
%
% RETURNS:
%   Phi - Reservoir pressure basis functions represented as a sparse
%         matrix of size G.cells.num-by-nchc.

   lam  = accumarray(p, mob .* G.cells.volumes) ./ ...
          accumarray(p,        G.cells.volumes);
   dlam = CS.mobility(CS.activeCellFaces)       ./ ...
          lam(CG.cellFaces(CS.activeCellFaces, 1));

   nc  = numel(CS.activeCellFaces);
   Phi = CS.basisP(:, CS.activeCellFaces) * spdiags(dlam, 0, nc, nc);
end

    
function s = id(s)
   s = ['solveIncompFlowMS:', s];
end

%--------------------------------------------------------------------------

function solver = pick_solver(CG, S, CS, dF, dC, opt)
   regul = ~any(dF);  % Set zero level if no prescribed pressure values.

   switch lower(opt.Solver),
      case 'hybrid',
         if ~any(strcmp(S.type, {'hybrid', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                  ['Solver mode ''hybrid'' is not consistent with ', ...
                   'linear system type ''%s''.'], S.type);
         end

         solver = @solve_hybrid;
      case 'mixed',
         if ~any(strcmp(S.type, {'mixed', 'tpfa', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                  ['Solver mode ''mixed'' is not consistent with ', ...
                   'linear system type ''%s''.'], S.type);
         end

         solver = @(A,b) solve_mixed(A, b, @mixedSymm);
      case 'tpfa',
         if ~any(strcmp(S.type, {'mixed', 'tpfa', 'comp_hybrid'})),
            error(id('SolverMode:Inconsistent'), ...
                  ['Solver mode ''tpfa'' is not consistent with ', ...
                   'linear system type ''%s''.'], S.type);
         end

         solver = @(A,b) solve_mixed(A, b, @tpfSymm);
      otherwise,
         error(id('SolverMode:NotSupported'), ...
               'Solver mode ''%s'' is not supported.', opt.solver);
   end
   %-----------------------------------------------------------------------
   % Enter prescribed pressures into final solution vector ----------------
   %
   lam     = zeros(size(dF));
   lam(dF) = dC;

   function x = solve_hybrid(A, b)
      A{3} = A{3}(:,~dF);
      b{3} = b{3}(  ~dF);
      [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                     'regularize', regul, ...
                                     'LinSolve', opt.LinSolve);
      lam(~dF) = x{3};
      x{3}     = lam;
   end

   function x = solve_mixed(A, b, solver)
      actF = CS.activeFaces;
      actB = CS.activeCellFaces;

      nF = neumann_faces(CG, dF, actF, opt);
      
      orientation_c = 2*(CG.cellFaces(:,1) == ...
                      CG.faces.neighbors(CG.cellFaces(:,2),1)) - 1;
      orientation = orientation_c(actB);
      Do = oriented_mapping(CS, orientation, opt);
        
      A{3} = A{3}(:,nF);
      b{3} = b{3}(  nF);
      
      % must remove 'MixedB' = 'true' when basis is on cellFaces
      %
      [x{1:4}] = solver(A{:}, b{:}, Do, 'Regularize', regul, ...
                        'LinSolve', opt.LinSolve,'MixedB', true); 
   
      % put contact pressure - NB: only defined on boundary
      %
      lam(nF) = x{3};
      x{3}    = lam;
   end
end

%--------------------------------------------------------------------------

function [xr, xw] = pack_solution(xr, xw, G, CG, p, S, CS, flux, pres, ... 
                                  lam, Phi, Psi, fluid, opt)
   W = opt.wells;   
   
   Lt = fluid.Lt(xr);
   [BW, Psi_w, Phi_w, CW, DW, fW, hW] = ...
      unpackWellSystemComponentsMS(W, G, p, Lt);
   %%
   [nsub, sub] = get_subfaces(G, CG, CS);
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;
   
   %-----------------------------------------------------------------------
   %% Package reservoir (coarse) solution structure -----------------------
   %
   
   switch lower(opt.Solver),
      case 'hybrid',

         xr.blockPressure   = pres;
         xr.blockFlux       = flux(1 : numel(actB));
         xr.cellPressure(:) = pres(p)                  + ...
            Phi   * xr.blockFlux + ...
            Phi_w * (-flux(numel(actB) + 1 : end));

         xr.facePressure(:) = accumarray(G.cellFaces(:,1), ...
            S.C * xr.cellPressure) ./ ...
            accumarray(G.cellFaces(:,1), 1);

         xr.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);

         xr.cellFlux(:)     = Psi * flux;
         xr.faceFlux(:)     = cellFlux2faceFlux(G, xr.cellFlux);

         fluxW   = flux(numel(actB) + 1 : end);
      case 'mixed'

         orientation_c = 2*(CG.cellFaces(:,1) == ...
                         CG.faces.neighbors(CG.cellFaces(:,2),1)) - 1;
         orientation = orientation_c(actB);

         xr.blockPressure   = pres;
         % convert to from coarse faceFlux to cellFaceFlux..
         %
         resFlux =  flux(1 : numel(actF));
         xr.blockFlux = orientation.* resFlux(CG.cellFaces(actB,2));

         xr.cellPressure(:) = pres(p); %                  + ...
         xr.facePressure(:) = accumarray(G.cellFaces(:,1), ...
                              S.C * xr.cellPressure) ./ ...
                              accumarray(G.cellFaces(:,1), 1);

         xr.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);

         xr.cellFlux(:)    = Psi * flux;
         xr.faceFlux(:)    = cellFlux2faceFlux(G, xr.cellFlux);

         fluxW   = flux(numel(actF) + 1 : end);
   end
   
   %-----------------------------------------------------------------------
   %% Package well solution structure -------------------------------------
   %
   %fluxW   = flux(numel(actF) + 1 : end);
   lamW    = lam (numel(actF) + 1 : end);
  
   %nb: fungerere bare for mixed!
   xw = packageWellSol(fluxW, lamW, fW, hW);

   %-----------------------------------------------------------------------
   %% Recover fine-grid fluxes in wells. ----------------------------------
   %
   i = 0;
   for k = 1 : numel(xw),
      rates = W(k).CS.rates;

      xw(k).flux = - full(rates * fluxW(i + 1 : i + size(rates,2)));

      i = i + size(rates,2);   
   end   
end

%--------------------------------------------------------------------------

function nF = neumann_faces(CG, dF, actF, opt)
% Determine the 'Neumann faces' (nF) (i.e., the faces for which
% face/contact pressures will be determined in the 'mixed' and 'TPFA'
% cases) of the model (g) according to the following rules:
%
%   - No internal faces are counted amongst the 'Neumann faces'.
%   - An external face is a 'Neumann face' unless there is a prescribed
%     pressure value (i.e., a Dirichlet condition) associated to the face.
%   - Wells controlled by rate-constraints are treated as Neumann contacts.
%
   
   nF = false(numel(actF),1); % No int.  
   % All external active faces ... :
   nF(any(CG.faces.neighbors(actF,:) == 0, 2)) = true;  
   nF(dF(1:numel(actF))) = false; % ... except Dirichlet cond.

   if ~isempty(opt.wells),
      % Additionally include all rate-constrained wells.
      nF = [nF; strcmp('rate', { opt.wells.type } .')];
   end
end

%--------------------------------------------------------------------------

function Do = oriented_mapping(CS, orient, opt)
% Do is the matrix mapping face-fluxes to half-face-fluxes and is used as
% input to mixedSymm and tpfSymm

   if ~isempty(opt.wells),
      ws = [ opt.wells.CS ];
      dw = diag(vertcat(ws.D));
      ow = ones([size(dw,1), 1]);
   else
      dw = [];
      ow = [];
   end
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;
   nf = numel(orient) + numel(ow); 
   
   Do = spdiags([orient; ow], 0, nf, nf) * blkdiag(CS.D(actB, actF), dw);
end
