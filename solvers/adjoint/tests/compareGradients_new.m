%compareGradients - compare gradients computed numerically and with adjoint

initSimpleModel

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                              
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
grad   = computeGradient(W, adjRes, schedule, controls);  

numGrad = computeNumericalGradient(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction)
adjGrad = cell2mat(grad)


figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

