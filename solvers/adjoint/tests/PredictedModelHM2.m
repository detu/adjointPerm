% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM2
% a 10 x 1 reservoir for history matching

nx = 10; ny = 1; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
m     = ones(nx,ny);
% Kreal = 2*ones(nx,ny);
Kreal = rand(10,1);
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;

perm = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

modelHMSmall2


% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                         
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
                
grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction);

% Compare with Finite Difference Grad
numGrad = computeNumericalGradientDK2(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction);
                            
adjGrad = cell2mat(grad);
adjGrad = adjGrad';

figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

% Start Optimization
                          
return;                              
                              
 