% compareGradients - compare gradients computed numerically and with adjoint
% for permeability field 
% Eka Suwartadi , 02-Nov-2009

initSimpleModel

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                              
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls);
grad   = grad';

numGrad = computeNumericalGradientDK(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction)
adjGrad = cell2mat(grad)
adjGrad = adjGrad';


figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

