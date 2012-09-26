% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% a 5 x 5 reservoir for history matching

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
m     = [ 0.7577    0.7060    0.8235    0.4387    0.4898; ...
          0.7431    0.0318    0.6948    0.3816    0.4456; ...
          0.3922    0.2769    0.3171    0.7655    0.6463; ...
          0.6555    0.0462    0.9502    0.7952    0.7094; ...
          0.1712    0.0971    0.0344    0.1869    0.7547 ];
      
Kreal = [ 1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1 ];

m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

% Choose objective function
objectiveFunction = str2func('FwMatchSmallReg');


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
                        

% -- OPTIMIZATION PART --
% initial permeability value
U   = param.m;
lb  = [];
ub  = [];

lbc = [];
ubc = [];

A   = []; 
b   = [];

Aeq = [];
beq = [];


auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param};
options = optimset('MaxIter', 20);

tic;
[U,fval,exitflag,output,lambda] = ktrlink(@(U)MisMatchWaterCut(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options,'knitro.opt');
toc;

fprintf('Elapsed time for optimizatin: %0.3f secs\n', toc);

% clf;
figure(1);
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
                              
 
