% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% a 4 x 4 reservoir for history matching

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% let K = one
m     = ones(nx,ny);
% m     = rand(5,5);
% m     = [ 1 2 3 4 5; ...
%           5 4 3 2 1; ...
%           1 2 3 4 5; ...
%           5 4 3 2 1; ...
%           1 2 3 4 5 ];
% m  = [0.1206 0.2518 0.9826 0.9063 0.0225; ... 
%       0.5895 0.2904 0.7302 0.8796 0.4252; ...
%       0.2261 0.6170 0.3438 0.8177 0.3127; ... 
%       0.3846 0.2652 0.5840 0.2607 0.1614; ...
%       0.5829 0.8243 0.1077 0.5943 0.1787];

% m  = [0.1 0.2 0.9 0.9 0.01; ... 
%       0.5 0.2 0.7 0.8 0.4; ...
%       0.2 0.6 0.3 0.8 0.3; ... 
%       0.3 0.2 0.5 0.2 0.1; ...
%       0.5 0.8 0.1 0.5 0.1];

Kreal = [ 1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1 ];
      
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

% save Kreal.mat Kreal m;
save Kreal.mat Kreal;
perm = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

% Choose objective function
objectiveFunction = str2func('FwMatchSmallReg');
%objectiveFunction = str2func('FwMatchSmall');

modelHMSmall3

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 4, 'TotalTime', 100*day, 'Verbose', verbose);
% schedule = initSchedule(W, 'TimeSteps', 30:30:30, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 4);
% parameters to be estimated
param.K = reshape(Kreal,G.cells.num,1);
param.m = repmat(m',numel(schedule),1);

% forward run
simRes = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, ...
                             'VerboseLevel', verboseLevel);
                         
% adjoint run
adjRes = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, param, 'VerboseLevel', verboseLevel);
                
grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);

% Compare with Finite Difference Grad
numGrad = computeNumericalGradientDK3(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction, param);
                            
adjGrad = cell2mat(grad);
adjGrad = 2*adjGrad';

adjGrad = repmat(adjGrad, 1, numel(schedule));

figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

% return;

% -- OPTIMIZATION PART --
% initial permeability value
U        = param.m;

% bound constraint
lb  = 1e-3*ones(100,1);
ub  = 1e2*ones(100,1);

opts = lbfgs_options('iprint', 1, 'maxits', 1000, 'factr', 1e5,'cb', @test_callback);

auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param};

tic;
[U,fx,exitflag,userdata] = lbfgs(@(U)MisMatchWaterCut(U,auxdata),U,opts);

toc;

fprintf('Elapsed time for optimizatin: %0.3f secs\n', toc);
                          
return;                              
                              
 
