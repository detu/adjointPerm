function numGrad = computeNumericalGradient(simRes, G, S, W, rock,     ...
                                            fluid, schedule, controls, ...
                                            objectiveFunction)
% compute numerical gradient

epsilon = 1e-5;

numControlWells = numel( controls.well );
numSteps        = numel( schedule);
obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
valInit   = obj.val;

uInit    = [controls.well.values]';
uInit    = uInit(:);
dimU     = length( uInit );

% scale epsilon
epsilon = epsilon*norm(uInit);

for k = 1:dimU
    e_k       = zeros(dimU, 1); e_k(k) = 1;
    uCur      = uInit + epsilon*e_k;
    controls  = updateControls(controls, uCur);
    schedule  = updateSchedule(controls, schedule);
    simRes    = runSchedule(simRes(1).resSol, G, S, W, rock, fluid, ...
                            schedule, 'VerboseLevel', 0);
    obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
    values(k) = obj.val;
end


numGrad = reshape((values' - valInit)/epsilon, numControlWells, numSteps);

function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end    