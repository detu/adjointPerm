function grad = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param)
% compute gradient for objective function with respect to permeability

DLtInv = @(sol)(-fluid.dkr(sol)*(1./fluid.mu)'...
                ./((1-sum(fluid.sr))*fluid.Lt(sol).^2));

numW  = numel(W);
numCF = size(S.B, 1);
obj   = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule, controls);
k     = numel(schedule);   % at final time step

l_v      = adjRes(k+1).resSol.cellFlux;
v        = simRes(k+1).resSol.cellFlux;
dim      = DLtInv(simRes(k+1).resSol);
cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C      = sparse(1:numel(cellNo), cellNo, 1);
grad{k}  = S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.DK*l_v;
grad{k}  = grad{k} + obj.partials(k+1).u';

for wellNr = 1:numW
    w       = W(wellNr);
    q       = simRes(k+1).wellSol(wellNr).flux; % q_w^{n+1}
    l_q     = adjRes(k+1).wellSol(wellNr).flux; % lam_q_w^{n+1}
    grad{k} = grad{k} - w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.DK*l_q) );
end


