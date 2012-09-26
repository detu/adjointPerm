function [Hv] = Hessvec(dv, u, usr_par)
% function [Hv] = Hessvec(G, S, W, rock, fluid, simRes, adjRes, schedule, controls, objectiveFunction, varargin)
% hessVec -- Compute Hessian times vector

 verboseLevel = 0;
[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction] = deal(usr_par{:});

m         = u;
param.m   = m;

% FORWARD SOLVE
simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);


% ADJOINT SOLVE
adjRes    = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction);

linSimRes = runLinearizedPerm1(simRes, G, S, W, rock, fluid, schedule, dv,  'VerboseLevel', 0);

secAdjRes = runSecondOrderAdjointPerm(simRes, adjRes, linSimRes, G, S, W, rock, fluid, ...
                                  schedule, controls, objectiveFunction, param, modelFunction, 'VerboseLevel', verboseLevel);

Hv        = compHvPerm(param, G, S, W, rock, simRes, secAdjRes, schedule, controls, fluid, objectiveFunction, dv); 

return
                              
%--------------------------------------------------------------------------

function Hv = compHvPerm(param, G, S, W, rock, simRes, secAdjRes, schedule, controls, fluid, objectiveFunction, dv)
% contribution from the obj.func:  partials(step).u2
% from the model: l_v and l_qw

DLtInv = @(sol)(-fluid.dkr(sol)*(1./fluid.mu)'...
                ./((1-sum(fluid.sr))*fluid.Lt(sol).^2));
numW   = numel(W);
numCF  = size(S.B, 1);
obj    = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule, controls);

k      = 1;
l_v    = secAdjRes(k+1).resSol.cellFlux;
v      = simRes(k+1).resSol.cellFlux;
dim    = DLtInv(simRes(k+1).resSol);
cellNo = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C    = sparse(1:numel(cellNo), cellNo, 1);
Hv     = S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.DK*l_v;
Hv     = Hv + (obj.partials(k+1).u2 .* dv);
% Hv     = Hv + (obj.partials(k+1).u2 );

for wellNr = 1:numW
    w       = W(wellNr);
    q       = simRes(k+1).wellSol(wellNr).flux; % q_w^{n+1}
    l_q     = secAdjRes(k+1).wellSol(wellNr).flux; % lam_q_w^{n+1}
    Hv      = Hv - w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.DK*l_q) );
end

return
                             