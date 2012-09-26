% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% a 4 x 4 reservoir for history matching

nx = 4; ny = 4; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% let K = one
% m     = ones(nx,ny);
m     = [ 2 2 3 3; ...
          2 2 3 3; ...
          4 4 1 1; ...
          4 4 1 1 ];
% m     = rand(4,4);
Kreal = [ 3 7 1 1; ...
          4 8 1 1; ...
          5 9 1 1; ...
          6 1 1 1 ];
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal m;

% load permSmallInit;
% perm = reshape(permSmallInit', [1,G.cells.num]);
perm = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
% rock.perm = perm'*100*milli*darcy;
% rock.perm = perm';

modelHMSmall

% nx = 15; ny = 15; nz = 1;
% G  = cartGrid([nx ny nz]);
% G  = computeGeometry(G);
% 
% load predPerm;
% perm = reshape(perm', [1,G.cells.num]);
% rock.perm = perm'*100*milli*darcy;
% 
% modelHM


% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                         
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
                
grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction);

% Compare with Finite Difference Grad
numGrad = computeNumericalGradientDK(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction);
                            
adjGrad = cell2mat(grad);
% adjGrad = -adjGrad';
% adjGrad = vertcat(adjGrad(:,1),adjGrad(:,2));
adjGrad = adjGrad';

figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

% Start Optimization
                          
return;                              
                              
 