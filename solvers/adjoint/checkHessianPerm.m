function [hessNum, hess] = checkHessianPerm()
% check Hessian obtained by second order adjoint code by comparing to 
% numerical Hessian obtained by perturbing gradient (by adjoint) at 
% current u by epsilon*e_k , k = 1,...,dimU.  

clear;
clc;
addpath ../adjoint

% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% let K = one
% m     = ones(nx,ny);
m     = rand(5,5);

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
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

% Choose objective function
objectiveFunction = str2func('FwMatchSmallReg');
% objectiveFunction = str2func('FwMatchSmall');

modelHMSmall3

% totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 100*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);

% parameters to be estimated
param.K = reshape(Kreal,G.cells.num,1);
param.m = repmat(m',numel(schedule),1);

% % forward run
% simRes = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, ...
%                              'VerboseLevel', verboseLevel);
%                          
% % adjoint run
% adjRes = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, ...
%                     objectiveFunction, param, 'VerboseLevel', verboseLevel);
%                 
% grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
% 
% % Compare with Finite Difference Grad
% numGrad = computeNumericalGradientDK3(simRes, G, S, W, rock, fluid, ...
%                                    schedule, controls, objectiveFunction, param);
%                             
% adjGrad = cell2mat(grad);
% adjGrad = adjGrad';
% 
% adjGrad = repmat(adjGrad, 1, numel(schedule));

% figure(1); hold on
% for k = 1 : size(numGrad, 1);
%     plot(numGrad(k,:), '-ob');
%     plot(adjGrad(k,:), '-xr');
% end
% legend('Numerical', 'Adjoint')

% return;

hn  = computeNumericalHessianPerm(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param);
h   = computeHessianPerm(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param);


disp(['Relative norm of differerence in Hessians: ' num2str(norm(hn-h)/norm(hn))]);
disp(['Relative norm of non-symmetric part: ' num2str(norm(h-h')/norm(h))]);

hessNum = hn;
hess    = h;