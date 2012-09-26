% animation of setupCase1
% comparing BFGS result against TN result

clear all

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

% % whether or not to show output
% verbose = false;
% verboseLevel = 0;
% 
% [G, S, rock, W, fluid] = initModel(2, 1);
% 
% resSolInit    = initResSol(G, 0.0);
% resSolInit.sw = 0.2*ones(G.cells.num, 1);
% 
% 
% % Choose objective function
% objectiveFunction = str2func('simpleNPVP');
% 
% % Hessian-vector product function
% HessianFunction = str2func('HessMultFcn');
% 
% totRate = 2;
% 
% totVol   = sum(G.cells.volumes.*rock.poros);
% schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', totVol, 'Verbose', verbose);
% 
% 
% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                        'Verbose', verbose,...
%                        'NumControlSteps', 5);

% load BFGS\TN result
% load BFGSResults;
load TNResults;

% Objective function NPV
% controls  = updateControls(controls, Ubfgs);
controls  = updateControls(controls, Utn);
schedule  = updateSchedule(controls, schedule);

% FORWARD SOLVE
simRes    = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

% Compute NPV
[obj,tid] = objectiveFunction(G, S, W, rock, fluid, simRes);

% AVI setting
num_frames_per_second = 2;
clear mex;
% aviobj = avifile('BFGSres.avi','fps', .6 ,'compression','indeo5');
aviobj = avifile('TNres.avi','fps', .6 ,'compression','indeo5');
aviobj.Quality = 100;

% plot oil saturation
% clf;
% fig = figure;
% set(fig,'DoubleBuffer','on');
% mov = avifile('BFGSres.avi');
% h = figure('visible','off');
h = figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(2,2,1);
plotCellData(G, log10(rock.perm));
axis off; 
plotWell(G, W, 'radius', .1*10, 'height', 2*10, 'color', 'k');
view([0 72]);
title('Permeability Field');

% use subplot
% 1) Permeability, 2) Oil Saturation Movement
% 3) Control Input,4) NPV per time

% prepare U, control input, to be displayed
Uinit = [4;-.25*4*ones(4,1)];  % read manually from input file
% Uall  = vertcat(Uinit,Ubfgs);
Uall  = vertcat(Uinit,Utn);
Ucon  = reshape( Uall, numel(schedule), numel(Uall)/numel(schedule)  ); 

% prepare NPV
val  = vertcat(0,tid);
times= [schedule.timeInterval]';
time = unique(times);

for k = 1: numel(simRes)
    
    subplot(2,2,2,'replace');
    
    if k == 1
        hs = plotCellData(G, 1-simRes(k).resSol.sw,(1 : G.cells.num),'FaceColor','r');
    else
        hs = plotCellData(G, 1-simRes(k).resSol.sw);
    end
    axis off; 
    plotWell(G, W, 'radius', .1*10, 'height', 2*10, 'color', 'k');
    view([0 72]);
    title('Oil Saturation');
    
    subplot(2,2,3,'replace');
    barh(Ucon(1,k),'facecolor','b');
    hold on;
    barh([0;-Ucon(2:5,k)],'facecolor','r');
    legend('Injector','Producer')
    axis off;
    title('Well Rates/Control Inputs');
    
    subplot(2,2,4,'replace');
    plot(time(1:k), val(1:k),'-ob'); grid
    set(gca,'Box','off');
    title('Obj.function: NPV');
    
%     frame = getframe ( gca );
    set(gca,'Zlim',[-20 20]); 
    frame = getframe ( h );
    aviobj = addframe ( aviobj, frame );
end

aviobj = close ( aviobj );
return;