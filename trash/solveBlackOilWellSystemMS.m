function [resSol, wellSol] = solveBlackOilWellSystemMS(resSol, wellSol, G, CG, rock, ...
                                                       S, CS, W, ...
                                                       PVTTAB, p0, dt, varargin)

fluid = pvt(PVTTAB, rock, resSol.cellPressure, S.constants.gravity);
CS = updateBlackOilMimeticSystemMS(G, CG, S, CS, rock, fluid, p0, resSol, dt,...
                                   varargin{:});
W  = updateBlackOilWellSystemMS   (G, CG, S, CS, W,    fluid, wellSol, ...
                                   varargin{:});


[BW, basisW, CW, DW, fW, hW] = unpackWellSystemComponentsMS(W, fluid.Lt);

%% Build full system
activeB = CS.activeCellFaces;
activeF = CS.activeFaces;
C = vertcat(CS.C(activeB,:),        CW{:});
D = blkdiag(CS.D(activeB, activeF), DW{:});

WCS = [W.CS];
BIV = vertcat(CS.BIV, WCS.BIV);
P = CS.P;

%  BC is really missing...
f = vertcat(CS.RHS.f_bc(activeB),      fW{:});
g = CS.RHS.g_src + CS.RHS.g_comp + CS.RHS.volume_discrepancy;
h = vertcat(CS.RHS.h_bc(activeF),      hW{:});

%% Finally Bms:
fullBasis = [CS.basis(:, activeB), -horzcat(basisW{:})];
B         = fullBasis' * (spdiags(S.C * (1./fluid.Lt), 0, S.sizeB(1), S.sizeB(2))*S.BI) * fullBasis;

wellRange = length(CS.activeCellFaces) + 1 : size(B,1);
B(wellRange, wellRange) = B(wellRange, wellRange) + blkdiag(BW{:});

%% Invert:
BI = sparse(size(B,1), size(B,1));
for k = 1 : size(C,2),
   range = find(C(:,k));
   BI(range, range) = inv( B(range, range) );
end

%% Solve non-symmetric Schur complement system.
dirFR = CS.RHS.dirichletFaces(activeF);
if ~isempty(W),
   dirFW = strcmp({ W.type } .', 'bhp');
else
   dirFW = logical([]);
end
dirF  = [dirFR; dirFW];

[flux, press, lam] = schurComplement(BI, C, D, BIV, P, f, g, h, dirF);


%% Package new reservoir solution
resSol.blockPressure   = press;
resSol.blockFlux       = flux;
resSol.cellPressure(:) = coarse2fine(press, CG);
resSol.cellFlux        = S.BI * (fullBasis * flux);
resSol.faceFlux        = cellFlux2faceFlux(G, resSol.cellFlux);

lamW(~dirFW) = lam(numel(activeF) + 1 : end);
if ~isempty(dirFW),
   lamW(dirFW) = [ W(dirFW).val ] .';
end
wellSol = packageWellSol(flux(numel(activeB) + 1 : end), ...
                         lamW, fW, hW);

% Recover fine-grid fluxes in wells
for w = 1 : numel(wellSol),
   wellSol(w).flux = S.C(:, W(w).cells)' * resSol.cellFlux;
end
