function [resSol, wellSol] = solveMixedWellSystemMS(G, S, CG, CS, rock, W, varargin)
% solveMixedWellSystemMS -- Solves coarse mixed system given B, C, D and RHS
%
% SYNOPSIS:
%   [resSol, wellSol] = solveMixedWellSystemMS(G, S, CG, CS, rock, W)
%
% PARAMETERS:
%   G       - SAMSIM grid data structure.
%
%   S       - System structure as defined by assembleMimeticSystem &c.
%
%   CG      - Coarse grid as defined by generateCoarseGrid.
%
%   CS      - Coarse system structure as defined by generateCoarseSystem.
%
%   rock    - Rock data structure.
%
%   W       - Well system structure as defined by assembleWellSystem and
%             generateCoarseWellSystem.
%
% RETURNS:
%   resSol  - Reservoir solution structure.
%
%   wellSol - Well solution structure.
%
% SEE ALSO:
%   solveMixedSystem, solveMixedWellSystem, inflow.

% $Id: solveMixedWellSystemMS.m 1408 2009-02-23 10:36:25Z ilig $

opt = struct('TotalMobility', []);
opt = merge_options(opt, varargin{:});
totMobFunc = opt.TotalMobility;

% Assemble well sub matrices
numWells = length(W);
if ~isempty(totMobFunc),
   if isa(totMobFunc, 'function_handle'),
      mobfunc = totMobFunc;
   elseif ischar(totMobFunc),
      mobfunc = str2func(totMobFunc);
   else
      error(id('TotalMobility:Unrecognized'), ...
            'Total mobility function of this type is not supported');
   end
   mob = mobfunc(rock.saturation);
else
   mob = ones([G.cells.num, 1]);
end
[BW, basisW, CW, DW, fW, hW] = unpackWellSystemComponentsMS(W, mob);

% Build full system
activeCF = CS.activeCellFaces;
activeF  = CS.activeFaces;
numACF   = numel(activeCF);
numAF    = numel(activeF);

C    = vertcat(CS.C(activeCF,:),        CW{:});
D    = blkdiag(CS.D(activeCF, activeF), DW{:});
DHat = blkdiag(CS.D(activeCF, activeF), diag(vertcat(DW{:})));
f    = vertcat(CS.RHS.f_bc(activeCF),   fW{:});
g    = CS.RHS.g_src;
h    = vertcat(CS.RHS.h_bc(activeF),       hW{:});

% Finally Bms:
orientation = G.cellFaces(:,2); 
numcf       = numel(orientation);
Do          = spdiags(orientation, 0, numcf, numcf) * S.D;

fullBasis = Do * [CS.basis(:, activeF), -horzcat(basisW{:})];
B         = fullBasis' * ( spdiags(S.C*(1./mob), 0, S.sizeB(1), S.sizeB(2))*S.B ) * fullBasis;

wellRange = numel(CS.activeFaces) + 1 : size(B,1);
B(wellRange, wellRange) = B(wellRange, wellRange) + blkdiag(BW{:});

orientation_c = 2*(CG.cellFaces(:,1) == CG.faces.neighbors(CG.cellFaces(:,2),1)) - 1;
orientation = [orientation_c(activeCF); ...
               ones(size(DHat,1)-nnz(activeCF), 1)];

neuFR = CS.RHS.neumannFaces(activeF);
neuFW = strcmp({ W.type } .', 'rate');
neuF  = [neuFR; neuFW];
              
%-----------------------------------------------------------------------
% Solve coarse mixed system --------------------------------------------
%
[fv, p, lam_n] = mixedSymm(B, C, D, f, g, h,  ...
                           orientation, neuF, ...
                           'MixedB', true, 'Face2CellFace', DHat);

%-----------------------------------------------------------------------
% Package solution in accessible form ----------------------------------
%
numNFR = nnz(CS.RHS.neumannFaces(activeF));
cellFlux = fullBasis * fv;
faceFlux = cellFlux2faceFlux(G, cellFlux);

resSol  = struct('cellPressure', p,                 ...
                 'facePressure', lam_n(1 : numNFR), ...
                 'faceFluxMS',   fv(1 : numAF),     ... 
                 'faceFlux',     faceFlux,          ...
                 'cellFlux',     cellFlux);
                 
             
fluxW        = fv(numAF + 1 : end);
lamW( neuFW) = lam_n(numNFR + 1 : end);
lamW(~neuFW) = [ W(~neuFW).val ] .';

wellSol = packageWellSol(fluxW, lamW, fW, hW);

% Recover fine-grid fluxes in wells
for w = 1 : numWells,
   wellSol(w).flux = S.C(:, W(w).cells)' * resSol.cellFlux;
end

%-----------------------------------------------------------------------

function s = id(s)
s = ['solveMixedWellSystemMS:', s];
