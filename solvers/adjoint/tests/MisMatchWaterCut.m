function [f,g] = MisMatchWaterCut(U,auxdata)
% note: need to save results from adjoint's run
% in order to speed-up CPU time !

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction] = deal(auxdata{:});


% arrange U/perm
% param.m = repmat(U,numel(schedule),1);
param.m = U;

% FORWARD SOLVE
simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);

% COMPUTE NPV
obj = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule);
f   = obj.val;

if nargout > 1
    
    % ADJOINT SOLVE
    adjRes = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction);
    save simAdj.mat simRes adjRes;
    
    % COMPUTE GRADIENT
    grad = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
    g    = cell2mat( {grad{:}}');
    g    = g';
    g    = repmat(g, 1, numel(schedule));

end
end
    
