function [f,g] = NPVfunction(U, auxdata)
% NPV/Objective function 

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction] = deal(auxdata{:});

% Objective function NPV
controls  = updateControls(controls, U);
schedule  = updateSchedule(controls, schedule);

% FORWARD SOLVE
simRes    = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

% COMPUTE NPV
%[obj,heval] = objectiveFunction(G, S, W, rock, fluid, simRes);
obj         = objectiveFunction(G, S, W, rock, fluid, simRes);
f           = -obj.val;
%h           = heval;

%if nargout > 2
if nargout > 1
    
    % ADJOINT SOLVE
    adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls,objectiveFunction);
    
    %save simAdj.mat simRes adjRes;
    
    % COMPUTE GRADIENT
    grad = computeGradient(W, adjRes, schedule, controls);
    g    = -cell2mat( {grad{:}}');

end
    

function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end