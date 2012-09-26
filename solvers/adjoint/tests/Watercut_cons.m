function [c ceq gradc gradceq]=Watercut_cons(U,auxcons)
% Water cut constraint and its derivative

rho = 1e-4;

% no nonlinear equality constraint
ceq     = [];
% gradceq = [];

[G, S, W, rock, fluid, resSolInit, schedule, controls, constrainFunction, param, modelFunction] = deal(auxcons{:});

param.m = U;

% FORWARD SOLVE
simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);

% COMPUTE CONSTRAINT, remember MATLAB form: c(x) <= 0
obj = constrainFunction(param, G, S, W, rock, fluid, simRes, schedule);
% c   = rho - obj.val;
c   =  obj.val - rho;

if nargout > 2
    adjRes = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, constrainFunction, param, modelFunction);
    gradc  = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, constrainFunction, param);
    gradc  = cell2mat( {gradc{:}}');
    gradceq= [];
end

return;

