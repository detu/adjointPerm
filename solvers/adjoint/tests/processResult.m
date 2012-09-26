clear, clc;
warning off all;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

addpath ../adjoint
addpath ../optimization

% load MAT file from simulation result
% load numCase4Results;
load numCase4Hv;
% load numCase4Result;

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
simRes_refSmall = numCase4True(numCtrlIntv);
% optimization results
K         = Kreal ./ U;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
numControlInterval = 12;
[W, fluid, S, resSolInit, verbose, verboseLevel] = modelFunction(G, rock);
schedule = initSchedule(W, 'NumSteps', numCtrlIntv, 'TotalTime', 360*day, 'Verbose', verbose);
simRes   = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction,  ...
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
legend('Measurement', 'Estimated', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd1');

% Plot at Prd 2
subplot(2,2,2);
plot(interval(2:end),dObsAll(2,1:12), '-ob', interval(2:end), gmAll(2,1:12), '-xr');grid;
legend('Measurement', 'Estimated', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd2');

% Plot at Prd 3
subplot(2,2,3);
plot(interval(2:end),dObsAll(3,1:12), '-ob', interval(2:end), gmAll(3,1:12), '-xr');grid;
legend('Measurement', 'Estimated', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd3');

% Plot at Prd 4
subplot(2,2,4);
plot(interval(2:end),dObsAll(4,1:12), '-ob', interval(2:end), gmAll(4,1:12), '-xr');grid;
legend('Measurement', 'Estimated', 'Location', 'best')
xlabel('Time [day]');
ylabel('Water Cut');
title('Water Cut at Prd4');
