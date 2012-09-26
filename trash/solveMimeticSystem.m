function resSol = solveMimeticSystem(resSol, G, S, rock, fluid, varargin)
% solveMimeticSystem -- Solve mimetic pressure linear system
%                       for incompressible two-phase models without wells.
%
% SYNOPSIS:
%   resSol = solveMimeticSystem(resSol, G, S, rock, fluid)
%   resSol = solveMimeticSystem(resSol, G, S, rock, fluid, 'pn', pv, ...)
%
% PARAMETERS:
%   resSol  - Initialized reservoir solution structure containing valid
%             saturations and pressures.
%
%   G       - Grid structure as described in grid_structure.
%
%   S       - System structure.
%
%   rock    - Rock data structure with valid field 'perm'.
%
%   fluid   - Fluid object as defined by function initSimpleFluid.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are
%               - Verbose      -- Whether or not to time the computational
%                                 process and emit informational messages
%                                 during the computation.
%                                 Logical.  Default value: Verbose = FALSE.
%
%               - MatrixOutput -- Whether or not to return the coefficient
%                                 matrix used in computing the face
%                                 pressures.
%                                 Logical.  Default value: FALSE.
%
%               - LinSolve     -- Solver for the resulting linear systems.
%                                 Handle to a function supporting the
%                                 syntax
%                                        x = LinSolve(A,b)
%                                 for solving the system Ax=b of linear
%                                 equations.
%                                 Default value: LinSolve = @mldivide
%                                                (i.e., backslash).
%
%               - bc           -- Boundary condtion structure as defined by
%                                 function 'addBC'.  This structure
%                                 accounts for all external boundary
%                                 contributions to the reservoir flow.
%                                 Default value: bc = [] meaning all
%                                 external no-flow (homogeneous Neumann)
%                                 conditions.
%
%               - src          -- Explicit source contributions as defined
%                                 by function 'addSource'.
%                                 Default value: src = [] meaning no
%                                 explicit sources exist in the model.
%
% RETURNS:
%   resSol - Updated reservoir solution structure having new values for the
%            fields:
%              - cellPressure  : Pressure in all reservoir cells.
%              - facePressure  : Pressure on all reservoir faces.
%              - cellFlux      : Out flux on local cell faces.
%              - faceFlux      : Flux over global, indexed faces
%                                corresponding G.faces.neighbors.
%              - A             : System matrix.
%                                Only returned if 'MatrixOutput' is set.
%
% SEE ALSO:
%   assembleMimeticSystem, initResSol, initSimpleFluid, addBC, addSource.
%
% NOTE:
%   This function is for solving hybrid systems arising when there are no
%   wells present in the model.  See function 'solveWellSystem' for how to
%   solve models containing wells.

% $Id: solveMimeticSystem.m 1773 2009-03-19 12:47:11Z bska $

assert (~isempty(rock));
assert (~isempty(fluid));

% There should be no compressible terms
assert (~isfield(S.RHS, 'g_comp') || all(S.RHS.g_comp == 0));

opt = struct('Verbose',      false,     ...
             'MatrixOutput', false,     ...
             'LinSolve',     @mldivide, ...
             'bc', [], 'src', []);
opt = merge_options(opt, varargin{:});

%--------------------------------------------------------------------------
% Build hybrid system components ------------------------------------------
%
if opt.Verbose,
   fprintf('Computing component matrices ... ')
   tic
end

BI = spdiags(S.C * fluid.Lt(resSol), 0, ...
             S.sizeB(1), S.sizeB(2)) * S.BI;
C  = S.C;
D  = S.D;

[f, g, h, dF, dC] = mimeticSysRHS(G, S, fluid.omega(resSol), ...
                                  opt.bc, opt.src);

tocif(opt.Verbose)

%--------------------------------------------------------------------------
% Solve global hybrid system (symmetric) ----------------------------------
%
if opt.Verbose,
   fprintf('Solving pressure system ... ')
   tic
end

nres = 3;
if opt.MatrixOutput, nres = nres + 1; end

regul = ~any(dF);  % Set pressure zero level if no prescr. pressure values.
[res{1 : nres}] = schurComplementSymm(BI, C, D(:,~dF), f, g, h(~dF), ...
                                      'regularize', regul,           ...
                                      'LinSolve'  , opt.LinSolve);

tocif(opt.Verbose)

%--------------------------------------------------------------------------
% Package solution in accessible form -------------------------------------
%
resSol.cellPressure(:)   = res{2};
resSol.facePressure(~dF) = res{3};
resSol.facePressure( dF) = dC    ;
resSol.cellFlux(:)       = res{1};
resSol.faceFlux(:)       = cellFlux2faceFlux(G, res{1});

if opt.MatrixOutput, resSol.A  = res{4}; end
