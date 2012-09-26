% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% a 15 x 15 reservoir for history matching
% Control Interval larger than 1 !

clear, clc;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

addpath ../adjoint

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% let K = one
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
% objectiveFunction = str2func('regOnly');

modelHMSmall3

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 100*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
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
adjGrad = adjGrad';

adjGrad = repmat(adjGrad, 1, numel(schedule));
clf;
figure(1); hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint', 'Location', 'best')

adjGrad  = adjGrad';
numGrad  = numGrad';
compare  = [adjGrad numGrad];
numError = adjGrad - numGrad;
error    = sqrt(numError'*numError);
fprintf(1, 'Numerical error: %0.3f \n', error);

return;

% -- OPTIMIZATION PART --
% initial permeability value
U        = param.m;

% first attempt UNCONSTRAINED CASE, without bound constraint
% lb = 1e-3*ones(16,1);
% ub = 1e2*ones(16,1);
% lb  = 1e-3*ones(100,1);
% ub  = 1e2*ones(100,1);
lb  = [];
ub  = [];

lbc = [];
ubc = [];

A   = []; 
b   = [];

Aeq = [];
beq = [];

% options = optimset('Display', 'iter', 'MaxIter', 50, 'Hessian', 'lbfgs', 'DerivativeCheck', 'on',...
%                    'TolFun', 1e-50, 'Algorithm', 'active-set', 'GradObj', 'on', 'TolCon', 1e-70);

auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param};
options = optimset('MaxIter', 20);

tic;
[U,fval,exitflag,output,lambda] = ktrlink(@(U)MisMatchWaterCut(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options,'knitro.opt');
% [U,fval,exitflag,output,lambda] = fmincon(@(U)MisMatchWaterCut(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],[]);
toc;

fprintf('Elapsed time for optimizatin: %0.3f secs\n', toc);

% clf;
figure(2);
U    = reshape(U, [1,G.cells.num]);
Kest = Kreal ./ U;
Kest = reshape(Kest,nx,ny);
subplot(1,2,1), imagesc(Kest); axis off;
title('Estimated Permeability');

% figure(2);
load Kmodel.mat;
K = reshape(K,nx,ny);
subplot(1,2,2),imagesc(K); axis off;
title('True Permeability');

% % To be plotted: 1) permeability 2) obj.function match: i.e water cut
% 
% for i=1:numSteps
%     % plot
%     figure(i);
%     Uopt = reshape(U(:,1),nx,ny);
%     imagesc(Uopt); figure(gcf)
% end
                          
return;                              
                              
 
