% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% CHANGES:
% 1. increase the control interval !
% 2. make the true model as a function call, no longer *.mat file again

clear, clc;
warning off all;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

addpath ../adjoint
addpath ../optimization

% run the TRUE model
numCtrlIntv = 1;
observCaseTrue(numCtrlIntv);

nx = 25; ny = 1; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% m     = rand(25,1);
load minit.mat
% m     = 1e-5 + ( 1e2 - 1e-5).*rand(25,1);
Kreal = ones(25,1);

m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

global numControlInterval;
numControlInterval = 1;
modelFunction      = str2func('modelObs');
[W, fluid, S, resSolInit, verbose, verboseLevel] = modelFunction(G, rock);

% Choose objective function
objectiveFunction = str2func('PermMismatch');
constrainFunction = str2func('Jac_Fwn');
% objectiveFunction = str2func('Jac_Fwn');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 360*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:2), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% parameters to be estimated
param.K = Kreal';   
param.m = m';

% % forward run
% simRes = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, ...
%                              'VerboseLevel', verboseLevel);
%                                                   
% % adjoint run
% adjRes = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, ...
%                     objectiveFunction, param, modelFunction, 'VerboseLevel', verboseLevel);
%                 
% grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
% 
% % Compare with Finite Difference Grad
% numGrad = computeNumericalGradientDK3(simRes, G, S, W, rock, fluid, ...
%                                    schedule, objectiveFunction, param, resSolInit, modelFunction);
%                             
% adjGrad = cell2mat(grad);
% adjGrad = adjGrad';
% 
% % adjGrad = repmat(adjGrad, 1, numel(schedule));
% clf;
% figure(1); hold on
% for k = 1 : size(numGrad, 1);
%     plot(numGrad(k,:), '-ob');
%     plot(adjGrad(k,:), '-xr');
% end
% legend('Numerical', 'Adjoint', 'Location', 'best')
% 
% adjGrad  = adjGrad';
% numGrad  = numGrad';
% compare  = [adjGrad numGrad];
% numError = adjGrad - numGrad;
% error    = sqrt(numError'*numError);
% fprintf(1, 'Numerical error: %0.3f \n', error);
% 
% return;

% -- OPTIMIZATION PART --
% initial permeability value
U   = param.m;
lb  = 1e-5*ones(25,1);
ub  = 1e2*ones(25,1);

lbc = [];
ubc = [];

A   = []; 
b   = [];

Aeq = [];
beq = [];

auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction};
auxcons = {G, S, W, rock, fluid, resSolInit, schedule, controls, constrainFunction, param, modelFunction};
options = optimset('MaxIter', 0);

tic;
[U,fval,exitflag,output,lambda] = ktrlink(@(U)ambiguityValue(U,auxdata),U,A,b,Aeq,beq,lb,ub,@(U)Watercut_cons(U,auxcons),options,'knitro.opt');
toc;
% tic;
% [U,fval,exitflag,output,lambda] = ktrlink(@(U)ambiguityValue(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options,'knitro.opt');
% toc;

fprintf('Elapsed time for optimizatin: %0.3f secs\n', toc);

figure(1);
U    = reshape(U, [1,G.cells.num]);
Kest = Kreal ./ U;
Kest = reshape(Kest,nx,ny);
subplot(1,2,1), imagesc(Kest); axis off;
title('Estimated Permeability');

load Kmodel.mat;
K = reshape(K,nx,ny);
subplot(1,2,2),imagesc(K); axis off;
title('True Permeability');

% CHECK CONSTRAINT VALUE
param.m = U';
simRes  = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);
obj     = constrainFunction(param, G, S, W, rock, fluid, simRes, schedule);
c       = obj.val;

          
return;                       