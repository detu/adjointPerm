function [resSol, wellSol] = solveMixedWellSystem(resSol, G, S, W, ...
                                                  fluid, varargin)
% solveMixedWellSystem -- Solves system given fields B, C, D, and RHS.
%
% SYNOPSIS:
%   [resSol, wellSol] = solveWellSystem(resSol, G, S, W, fluid)
%
% PARAMETERS:
%   resSol - Initialized reservoir solution structure containing valid
%            saturations and cell pressures.
%
%   G      - Grid structure as described in grid_structure.
%
%   S      - System structure as defined by function assembleMimeticSystem.
%
%   W      - Well structure as defined by function addWell &c.
%
%   fluid  - Initialized fluid model structure as defined by function
%            initSimpleFluid.
%
% RETURNS:
%   resSol  - Reservoir solution structure having the following fields:
%               - cellPressure -- Pressure in all active reservoir cells.
%               - facePressure --
%               - cellFlux     -- Local cell faces out-flux.
%               - faceFlux     -- Flux over global, indexed faces
%                                 corresponding to G.faces.neighbors .
%
%   wellSol - Well solution structure having the following fields:
%               - flux         -- Well perforation flux.
%               - pressure     -- Well perforation pressure.
%
% SEE ALSO:
%   assembleMimeticSystem, assembleWellSystem, initResSol, initSimpleFluid.

% $Id: solveMixedWellSystem.m 1773 2009-03-19 12:47:11Z bska $

%--------------------------------------------------------------------------
% Extract well system contributions ---------------------------------------
%
Lt = fluid.Lt(resSol);
[BW, CW, DW, fW, hW] = unpackMixedWellSystemComponents(W, Lt);

%--------------------------------------------------------------------------
% Build hybrid system components ------------------------------------------
%
B = blkdiag(spdiags(S.C * (1 ./ Lt), ...
                    0, S.sizeB(1), S.sizeB(2)) * S.B, BW{:});
C = vertcat(S.C, CW{:});
D = blkdiag(S.D, DW{:});

DHat = blkdiag(S.D, diag(vertcat(DW{:})));

omega     = fluid.omega(resSol);
f_mimetic = S.RHS.f_bc + ...
            S.RHS.f_grav .* (S.C * omega) .* S.constants.gravity;

f = vertcat(f_mimetic,  fW{:});
g = S.RHS.g_src;
h = vertcat(S.RHS.h_bc, hW{:});

orientation  = [G.cellFaces(:,2); ...
                ones([size(B,1) - S.sizeB(1), 1])];

neuFW  = strcmp({ W.type } .', 'rate');
neuFR  = S.RHS.neumannFaces;
numNFR = nnz(neuFR);
neuF   = [neuFR; neuFW];

%--------------------------------------------------------------------------
% Solve global hybrid system (symmetric) ----------------------------------
%
[faceFlux, press, lam_n] = mixedSymm(B, C, D, f, g, h,  ...
                                     orientation, neuF, ...
                                     'Face2CellFace', DHat);

%--------------------------------------------------------------------------
% Package solution in accessible form -------------------------------------
%

resSol.cellPressure(:)     = press;
resSol.facePressure(neuFR) = lam_n(1 : numNFR);
resSol.cellFlux(:)         = faceFlux2cellFlux(G, faceFlux(1 : S.sizeD(2)));
resSol.faceFlux(:)         = faceFlux(1 : S.sizeD(2));

lamW( neuFW) = lam_n(numNFR + 1 : end);
lamW(~neuFW) = [ W(~neuFW).val ] .';

wellSol = packageWellSol(faceFlux(S.sizeD(2) + 1 : end), lamW, fW, hW);
