function hess = computeNumericalHessian(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction)
% compute numerical Hessian by comparing gradient (by adjoint) at current u
% and u+epsilon*e_k , k = 1,...,dimU 

epsilon = 1e-6;

%-------------- Initial gradient ------------------------------------------
simRes  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
adjRes  = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', 0);
grad    = computeGradient(W, adjRes, schedule, controls); 
gradVec = cell2mat( {grad{:}}');

%-------------- Gradients for perturbed controls --------------------------

u    = [controls.well.values]';
u    = u(:);
dimU = length(u);

pertGradVecs = zeros(dimU, dimU);

h = waitbar(0,'Computing gradients for perturbed controls ...');
for k = 1:dimU
    e_k       = zeros(dimU, 1); e_k(k) = 1;
    u_k       = u + epsilon*e_k;
    controls  = updateControls(controls, u_k);
    schedule  = updateSchedule(controls, schedule);
    simRes    = runSchedule(simRes(1).resSol, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
    adjRes  = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', 0);
    grad    = computeGradient(W, adjRes, schedule, controls); 
    pertGradVecs(:, k) = cell2mat( {grad{:}}');
    waitbar(k/dimU,h)
end
close(h)

hess = ( pertGradVecs - gradVec*ones(1,dimU) )/epsilon;


function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end    