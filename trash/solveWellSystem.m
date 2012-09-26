function [resSol, wellSol] = solveWellSystem(resSol, G, S, W, fluid, ...
                                             varargin)
% solveWellSystem -- Solves system given fields BI, C, D, and RHS.
%
% SYNOPSIS:
%   [resSol, wellSol] = solveWellSystem(resSol, G, S, W, fluid)
%   [resSol, wellSol] = solveWellSystem(resSol, G, S, W, fluid, ...
%                                       'pn1', pv1, ...)
%
% PARAMETERS:
%   resSol  - Initialized reservoir solution structure containing valid
%             saturations and cell pressures.
%
%   G       - Grid structure as described in grid_structure.
%
%   S       - System structure.
%
%   W       - Well structure, can be empty (W = struct([]))
%
%   fluid   - Fluid object as defined by function initSimpleFluid.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are
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
%   resSol  - Reservoir solution structure having the following fields:
%               - cellPressure -- Pressure in all active reservoir cells.
%               - facePressure --
%               - cellFlux     -- Out flux on local cell faces.
%               - faceFlux     -- Flux over global, indexed faces
%                                 corresponding to G.faces.neighbors .
%
%   wellSol - Well solution structure having the following fields:
%               - flux         -- Well perforation flux.
%               - pressure     -- Well perforation pressure.
%
% SEE ALSO:
%   assembleMimeticSystem, assembleWellSystem, initResSol, initSimpleFluid.

% $Id: solveWellSystem.m 1773 2009-03-19 12:47:11Z bska $

%--------------------------------------------------------------------------
% Extract well system contributions ---------------------------------------
%

opt = struct('bc', [], 'src', [], 'LinSolve', @mldivide);
opt = merge_options(opt, varargin{:});

Lt = fluid.Lt(resSol);
[BIW, BIVW, CW, DW, fW, hW] = unpackWellSystemComponents(W, Lt);

%--------------------------------------------------------------------------
% Build hybrid system components ------------------------------------------
%

BI = blkdiag(spdiags(S.C * Lt, 0, ...
                     S.sizeB(1), S.sizeB(2)) * S.BI, BIW{:});
C  = vertcat(S.C, CW{:});
D  = blkdiag(S.D, DW{:});

[f, g, h, dF, dC] = mimeticSysRHS(G, S, fluid.omega(resSol), ...
                                  opt.bc, opt.src);

f = vertcat(f, fW{:});
h = vertcat(h, hW{:});

%--------------------------------------------------------------------------
% Solve global hybrid system (symmetric) ----------------------------------

%--------------------------------------------------------------------------
% Eliminate known contact and well pressures from linsys ------------------
%
if ~isempty(W),
   dFw = strcmp({ W.type } .', 'bhp');
   dCw = [ W(dFw).val ].';
else
   dFw = logical([]);
   dCw = [];
end

dF = [dF; dFw];
D  = D(:,~dF);
h  = h(  ~dF);

%--------------------------------------------------------------------------
% Enter prescribed pressures into final solution vector -------------------
%
lam     = zeros(size(dF));
lam(dF) = [dC; dCw];

%--------------------------------------------------------------------------
% Solve actual system of linear equations for remaining p/v values --------
%
regul = ~any(dF);  % Set pressure zero level if no prescr. pressure values.
[flux, press, lam(~dF)] = schurComplementSymm(BI, C, D, f, g, h,   ...
                                              'regularize', regul, ...
                                              'LinSolve', opt.LinSolve);

%--------------------------------------------------------------------------
% Package solution in accessible form -------------------------------------
%
resSol.cellPressure(:) = press;
resSol.facePressure(:) = lam (1 : S.sizeD(2));
resSol.cellFlux(:)     = flux(1 : S.sizeB(1));
resSol.faceFlux(:)     = cellFlux2faceFlux(G, flux(1 : S.sizeB(1)));

lamW    = lam(S.sizeD(2) + 1 : end);
wellSol = packageWellSol(flux(S.sizeB(1) + 1 : end), lamW, fW, hW);
