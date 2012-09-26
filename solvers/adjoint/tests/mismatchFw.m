function [f,g] = mismatchFw(U,auxdata)
% note: need to save results from adjoint's run
% in order to speed-up CPU time !

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov, priorPerm] = deal(auxdata{:});

%prior permeability harus di-fix ! CONSTANT VALUE !!!
%priorPerm = param.m;
param.m   = U;

% FORWARD SOLVE
simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);

% UPDATE g_mEns HERE !!!
% collected water cut as measurement
resSol    = simRes(asimStep+1).resSol;
wellSol   = simRes(asimStep+1).wellSol;
[wellRates, rateSigns] = getRates(W, wellSol);
wellCells = vertcat( W.cells );
wellSats  = resSol.s( wellCells );

f_w_all   = fluid.fw(resSol);
f_w       = f_w_all(wellCells);
f_o       = 1 - f_w;
injInx    = (rateSigns > 0);
prodInx   = (rateSigns < 0);
g_mEns    = f_w( prodInx );

% COMPUTE NPV
obj = objectiveFunction(G, W, fluid, g_mEns, g_m, asimStep, cov, priorPerm, param);
f   = -obj.val;

if nargout > 1
    
    % ADJOINT SOLVE
    adData = {g_mEns, g_m, asimStep, cov, priorPerm};
    adjRes = runAdjointPermKF(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction, adData);
    save simAdj.mat simRes adjRes;
    
    % COMPUTE GRADIENT
    grad = computeGradientDKF(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, adData);
    g    = cell2mat( {grad{:}}');
    g    = -g';
    g    = repmat(g, 1, numel(schedule));

end
end
    
