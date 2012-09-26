% optimization_driver for Five water-spot case
% for checking derivative: first and Hessian-times-vector product

clear;
clc;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
addpath ../optimization
addpath ../adjoint

% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = rand(nx,ny);
% m  = 100*ones(nx,ny);

Kreal = [ 1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1 ];

% Kreal = rand(5,5);
      
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
perm = reshape(K, [1,G.cells.num]);
% rock.perm = perm'*100*milli*darcy;
rock.perm = perm'*100*milli;

% Choose objective function
% objectiveFunction = str2func('FwMatchSmallReg');
% objectiveFunction = str2func('FwMatchSmall');
objectiveFunction = str2func('regOnly');

modelHMSmall3  %update permeability

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 100*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);

% parameters to be estimated
param.K = reshape(Kreal,G.cells.num,1);
param.m = repmat(m',numel(schedule),1);

% % perform derivative checks
usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param};
Uinit   = param.m;
Uinit   = Uinit(:);
deriv_check( Uinit, 1, usr_par);

return

