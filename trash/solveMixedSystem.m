function resSol = solveMixedSystem(resSol, G, S, rock, fluid, varargin)
% solveMixedSystem -- Solve mixed pressure linear system
%                     for models without wells.
%
% SYNOPSIS:
%   resSol = solveMixedSystem(resSol, G, S, rock, fluid)
%   resSol = solveMixedSystem(resSol, G, S, rock, fluid, 'pn', pv, ...)
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
%              - facePressure  : Pressure on external Neumann faces.
%              - cellFlux      : Out flux on local cell faces.
%              - faceFlux      : Flux over global, index faces
%                                corresponding to G.faces.neighbors.
%              - A             : System matrix.
%                                Only returned if 'MatrixOutput' is set.
%
% SEE ALSO:
%   assembleMimeticSystem, initResSol, initSimpleFluid, solveMimeticSystem.
%
% NOTE:
%   This function is for solving mixed systems arising from a two-point
%   flux discretization of an incompressible two-phase pressure equation in
%   a model with no wells.  Use function 'solveMixedWellSystem' to solve
%   models containing wells, though 'solveMixedWellSystem' does not support
%   TPFA discretizations.

% $Id: solveMixedSystem.m 1773 2009-03-19 12:47:11Z bska $

assert (~isempty(rock));
assert (~isempty(fluid));

% There should be no compressible terms
assert (~isfield(S.RHS, 'g_comp') || all(S.RHS.g_comp == 0));

if ~strcmp(S.type, 'mixed')
   error(id('Type:NotMixed'), ...
         'System has to be of type ''mixed'' to use solveMixedSystem')
end

opt = struct('Verbose',      false, ...
             'MatrixOutput', false, ...
             'src', [], 'bc', []);
opt = merge_options(opt, varargin{:});

verbose      = opt.Verbose;
matrixOutput = opt.MatrixOutput;

%--------------------------------------------------------------------------
% Build hybrid system components ------------------------------------------
%

if verbose,
   fprintf('Computing component matrices ... ')
   tic
end

Lt     = fluid.Lt(resSol);
cellno = rldecode((1 : G.cells.num) .', double(G.cells.numFaces));
B      = spdiags(1 ./ Lt(cellno), 0, S.sizeB(1), S.sizeB(2)) * S.B;
C      = S.C;
D      = S.D;

[f, g, h, dF, dC] = mimeticSysRHS(G, S, fluid.omega(resSol), ...
                                  opt.bc, opt.src);

tocif(verbose)

%--------------------------------------------------------------------------
% Solve global mixed system (symmetric) -----------------------------------
%
if verbose,
   fprintf('Solving pressure system ... ')
   tic
end

nF                                 = false([G.faces.num, 1]);
nF(any(G.faces.neighbors == 0, 2)) = true;
nF(dF)                             = false;

nres = 3;
if matrixOutput, nres = nres + 1; end

regul = ~any(dF);  % Set pressure zero level if no prescr. pressure values.
nCF   = size(G.cellFaces,1);
Do    = spdiags(G.cellFaces(:,2), 0, nCF, nCF) * D;

[res{1 : nres}] = mixedSymm(B, C, D(:,nF), f, g, h(nF), Do, ...
                            'Regularize', regul);
tocif(verbose)

%--------------------------------------------------------------------------
% Package solution in accessible form -------------------------------------
%
resSol.cellPressure(:)  = res{2};
resSol.facePressure(nF) = res{3};
resSol.facePressure(dF) = dC;
resSol.faceFlux(:)      = res{1};
resSol.cellFlux(:)      = faceFlux2cellFlux(G, res{1});

if matrixOutput, resSol.A = res{4}; end

%--------------------------------------------------------------------------

function s = id(s)
s = ['solveMixedSystem:', s];
