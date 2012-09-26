function [resSol, wellSol] = solveWellSystemMS(resSol, G, CG, p, S, CS, ...
                                               W, fluid, varargin)
% solveWellSystemMS -- Solve coarse (multiscale) pressure system.
%
% SYNOPSIS:
%  [xr, xw] = solveWellSystemMS(xr, G, CG, p, S, CS, W, fluid)
%  [xr, xw] = solveWellSystemMS(xr, G, CG, p, S, CS, W, fluid, ...
%                               'pn1', pv1, ...)
%
% PARAMETERS:
%   xr       - Initialized reservoir solution structure containing valid
%              saturations and cell pressures.
%
%   G, CG, p - Grid, coarse grid, and cell-to-block partition vector,
%              respectively.
%
%   S, CS, W - Linear system structure on fine grid (S) and coarse grid
%              (CS) as defined by functions 'assembleMimeticSystem' and
%              'generateCoarseSystem', respectively in addition to well
%              linear system structure (W) as defined by functions
%              'addWell', 'assembleWellSystem', and
%              'assembleCoarseWellSystem'.  An empty well linear system
%              structure (i.e., W=struct([])) is supported.  The resulting
%              model is interpreted as not containing any well driving
%              forces.
%
%   fluid    - Fluid object as defined by function initSimpleFluid.
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.
%              The supported options are:
%                 - bc  -- Boundary condtion structure as defined by
%                          function 'addBC'.  This structure accounts for
%                          all external boundary contributions to the
%                          reservoir flow.
%                          Default value: bc = [] meaning all external
%                          no-flow (homogeneous Neumann) conditions.
%
%                 - src -- Explicit source contributions as defined by
%                          function 'addSource'.
%                          Default value: src = [] meaning no explicit
%                          sources exist in the model.
%
%                 - LinSolve --
%                          Solver for the resulting linear systems.
%                          Handle to a function supporting the syntax
%                                    x = LinSolve(A,b)
%                          for solving the system Ax=b of linear equations.
%                          Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   xr - Updated reservoir solution structure.
%   xw - Well solution structure.
%
% SEE ALSO:
%   solveMimeticSystem, solveBlackOilWellSystem.

% $Id: solveWellSystemMS.m 1408 2009-02-23 10:36:25Z ilig $

   opt = struct('src', [], 'bc', [], 'LinSolve', @mldivide);
   opt = merge_options(opt, varargin{:});

   [nsub, sub] = get_subfaces(G, CG, CS);
   [I, J]      = coarsening_ops(G, CG, CS, sub, nsub);

   %-----------------------------------------------------------------------
   % Assemble well sub matrices -------------------------------------------
   %
   Lt = fluid.Lt(resSol);
   [BW, Psi_w, Phi_w, CW, DW, fW, hW] = ...
      unpackWellSystemComponentsMS(W, G, p, Lt);

   omega = fluid.omega(resSol);

   %-----------------------------------------------------------------------
   %% Build full (block) system of linear equations -----------------------
   %
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;

   flx_bas = [CS.basis(:, actB), -Psi_w];
   Psi     = S.BI * flx_bas;
   Phi     = res_basisP(G, CG, CS, p, Lt);

   C = vertcat(CS.C(actB,:),    CW{:});
   D = blkdiag(CS.D(actB,actF), DW{:});

   [f, g, h, dFr, dCr] = sys_rhs(G, S, omega, Psi(:, 1 : numel(actB)), ...
                                 I, J, opt.src, opt.bc);

   % Add well contributions to linear system.
   f = vertcat(f, fW{:});
   h = vertcat(h, hW{:});

   %-----------------------------------------------------------------------
   %% Assemble (inverse) coarse mass matrix -------------------------------

   %-----------------------------------------------------------------------
   % Original mass matrix Bc = psi.' * (Bf/mob) * psi ---------------------
   %
   Lti = spdiags(S.C * (1 ./ Lt), 0, S.sizeB(1), S.sizeB(2));
   B   = (Lti * flx_bas).' * Psi;

   wr = numel(actB) + 1 : size(B,1);
   B(wr,wr) = B(wr,wr) + blkdiag(BW{:});

   %-----------------------------------------------------------------------
   % Construct BI = INV(B) which is needed when solving discrete system ---
   %
   BI = invert(B, C);

   %-----------------------------------------------------------------------
   %% Solve symmetric Schur complement system -----------------------------

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

   dF = [dFr; dFw];
   D  = D(:,~dF);
   h  = h(  ~dF);

   %-----------------------------------------------------------------------
   % Enter prescribed pressures into final solution vector ----------------
   %
   lam     = zeros(size(dF));
   lam(dF) = [dCr; dCw];

   %-----------------------------------------------------------------------
   % Solve actual system of linear equations for remaining p/v values -----
   %
   do_reg = sum(dF) == 0;  % Need to set pressure zero level?
   [flux, pres, lam(~dF)] = schurComplementSymm(BI, C, D, f, g, h,    ...
                                                'regularize', do_reg, ...
                                                'LinSolve', opt.LinSolve);

   %-----------------------------------------------------------------------
   %% Package reservoir (coarse) solution structure -----------------------
   %
   resSol.blockPressure   = pres;
   resSol.blockFlux       = flux(1 : numel(actB));
   resSol.cellPressure(:) = pres(p)                  + ...
                            Phi   * resSol.blockFlux + ...
                            Phi_w * (-flux(numel(actB) + 1 : end));

   resSol.facePressure(:) = accumarray(G.cellFaces(:,1), ...
                                       S.C * resSol.cellPressure) ./ ...
                            accumarray(G.cellFaces(:,1), 1);

   resSol.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);

   resSol.cellFlux(:)     = Psi * flux;
   resSol.faceFlux(:)     = cellFlux2faceFlux(G, resSol.cellFlux);

   %-----------------------------------------------------------------------
   %% Package well solution structure -------------------------------------
   %
   fluxW   = flux(numel(actB) + 1 : end);
   lamW    = lam (numel(actF) + 1 : end);
   wellSol = packageWellSol(flux(numel(actB) + 1 : end), lamW, fW, hW);

   % Recover fine-grid fluxes in wells.
   i = 0;
   for k = 1 : numel(wellSol),
      rates = W(k).CS.rates;

      wellSol(k).flux = - full(rates * fluxW(i + 1 : i + size(rates,2)));

      i = i + size(rates,2);
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [nsub, sub] = get_subfaces(g, cg, cs)
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

   [nsub, sub] = subFaces(g, cg);
   sub_ix      = cumsum([0; nsub]);
   sub         = sub(mcolon(sub_ix(cs.activeFaces) + 1, ...
                            sub_ix(cs.activeFaces + 1)));
   nsub        = nsub(cs.activeFaces);
end

%--------------------------------------------------------------------------

function [I, J] = coarsening_ops(g, cg, cs, sub, nsub)
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

   I = cg.cells.subCells;
   J = sparse(sub, rldecode(cs.activeFaces(:), nsub), ...
              1, g.faces.num, cg.faces.num);
   J = J(:, cs.activeFaces);
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = sys_rhs(g, s, omega, Psi, I, J, src, bc)
% sys_rhs -- Evaluate coarse system right hand side contributions.
%
% SYNOPSIS:
%   [f, g, h, dF, dC] = sys_rhs(G, S, omega, Psi, I, J, src, bc)
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
%   mimeticSysRHS.

   % Evaluate fine-scale system right hand side contributions.
   %
   [ff, gg, hh, dF, dC] = mimeticSysRHS(g, s, omega, bc, src);

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

function Phi = res_basisP(g, cg, cs, p, mob)
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
%         matrix of size g.cells.num-by-nchc.

   lam  = accumarray(p, mob .* g.cells.volumes) ./ ...
          accumarray(p,        g.cells.volumes);
   dlam = cs.mobility(cs.activeCellFaces)       ./ ...
          lam(cg.cellFaces(cs.activeCellFaces, 1));

   nc  = numel(cs.activeCellFaces);
   Phi = cs.basisP(:, cs.activeCellFaces) * spdiags(dlam, 0, nc, nc);
end
