%compareGradients - compare gradients computed numerically and with adjoint

initSimpleModel

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                         
obj      = objectiveFunction(G, S, W, rock, fluid, simRes);
Jinit    = obj.val;
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
grad   = computeGradient(W, adjRes, schedule, controls);  

adjGrad = cell2mat(grad);

JnormJscale = .2 * (adjGrad ./ norm(adjGrad,inf));
uInit       = [controls.well.values]';
uCur        = uInit + JnormJscale;
controls    = updateControls(controls, uCur);
schedule    = updateSchedule(controls, schedule);

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
obj      = objectiveFunction(G, S, W, rock, fluid, simRes);
Jpert    = obj.val;
                              
% test to ONE, psi should be ONE if correct !
psi = ( Jpert - Jinit ) ./ (JnormJscale' * adjGrad); 

return;



