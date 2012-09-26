% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% CHANGES:
% 1. increase the control interval !
% 2. make the true model as a function call, no longer *.mat file again

clear, clc;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

addpath ../adjoint

% run the TRUE model
numCtrlIntv = 1;
numCase2True(numCtrlIntv);

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
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
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 360*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% parameters to be estimated
param.K = Kreal';   
param.m = m';

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
%                                    schedule, objectiveFunction, param, resSolInit);
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

load Kmodel.mat;
K = reshape(K,nx,ny);
subplot(1,2,2),imagesc(K); axis off;
title('True Permeability');

% plot water cut, increase control interval !
% reference model / the truth
numCtrlIntv     = 12;
simRes_refSmall = numCase2True(numCtrlIntv);
% optimization results
K         = Kreal ./ U;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
modelHMSmall3
schedule = initSchedule(W, 'NumSteps', numCtrlIntv, 'TotalTime', 360*day, 'Verbose', verbose);
simRes   = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, ...
                             'VerboseLevel', verboseLevel);
numSteps = numel(simRes);
interval = [simRes.timeInterval]';

gmAll    = [];
dObsAll  = [];
for step = 2 : numSteps
    % model
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;
    
    % measurement/reference
    resSolref  = simRes_refSmall(step).resSol;
    wellSolref = simRes_refSmall(step).wellSol;
    
    
    % model
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells ); 

    
    % measurement/reference
    [wellRates_ref, rateSigns_ref] = getRates(W, wellSolref);

    % model
    f_w_all   = fluid.fw(resSol);
    f_w       = f_w_all(wellCells);
    f_o       = 1 - f_w;
    injInx    = (rateSigns > 0);
    prodInx   = (rateSigns < 0);
    
    % measurement/reference
    f_w_all_ref   = fluid.fw(resSolref);
    f_w_ref       = f_w_all_ref(wellCells);
    f_o_ref       = 1 - f_w_ref;
    injInx_ref    = (rateSigns_ref > 0);
    prodInx_ref   = (rateSigns_ref < 0);
    
    % model - liquid rate at producer wells
    g_m     = f_w( prodInx ) ;
    gmAll   = [gmAll g_m];
    
    % measurement - liquid rate at producer wells
    d_obs    = f_w_ref( prodInx_ref );
    dObsAll  = [dObsAll d_obs];
end

interval = unique(interval);
interval = convertTo(interval,day);
figure(2);
% Plot at Prd 1
subplot(2,2,1);
plot(interval(2:end),dObsAll(1,1:12), '-ob', interval(2:end), gmAll(1,1:12), '-xr');grid;
legend('Observation', 'Optimization', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd1');

% Plot at Prd 2
subplot(2,2,2);
plot(interval(2:end),dObsAll(2,1:12), '-ob', interval(2:end), gmAll(2,1:12), '-xr');grid;
legend('Observation', 'Optimization', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd2');

% Plot at Prd 3
subplot(2,2,3);
plot(interval(2:end),dObsAll(3,1:12), '-ob', interval(2:end), gmAll(3,1:12), '-xr');grid;
legend('Observation', 'Optimization', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd3');

% Plot at Prd 4
subplot(2,2,4);
plot(interval(2:end),dObsAll(4,1:12), '-ob', interval(2:end), gmAll(4,1:12), '-xr');grid;
legend('Observation', 'Optimization', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd4');
                          
return;                              
                              
 
