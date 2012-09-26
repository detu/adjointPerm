% Using optimal well-rates at permCase1, perform measurement at producer
% wells
clear, clc;
warning off all;
opengl neverselect;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
% whether or not to show output
verbose = false;
verboseLevel = 0;

% load 101 realizations, unit is in milli darcy
load permeability;

% reservoir setting
nx = 45; ny = 45; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

% rock properties
m     = reshape(m, [1,G.cells.num]);
Kreal = permeability(:,101);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
save Kmodel.mat K;

perm      = reshape(K, [1,G.cells.num]);
perm      = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

% fluid properties
fluid         = initCoreyFluid('mu', [1, 5], 'sr', [0.1 0.1]);
resSolInit    = initResSol(G, 0.0);
resSolInit.s  = 0.2*ones(G.cells.num, 1);
S             = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');

% Introduce wells
radius  = .1;
totRate = 1/day;
W = [];
W = addWell(G, rock,W, 1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', 'I1');
W = addWell(G, rock,W, 45     , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', 'P1');
for i=1:5
    W = addWell(G, rock,W, (9*i - 1)*nx +1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', ['I',num2str(i+1)]);
end

for i=1:5
    W = addWell(G, rock,W, i*9*nx      , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', ['P',num2str(i+1)]);
end

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol     = sum(G.cells.volumes.*rock.poro);
schedule   = initSchedule(W, 'TimeSteps', 120*day:120*day:2400*day, 'Verbose', verbose);

injMaxMin  = repmat( [0.001 100]*meter^3/day, 6, 1);
prodMaxMin = repmat( [-100 -0.001]*meter^3/day, 6, 1);

controls   = initControls(schedule, 'ControllableWells', (1:12), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);
% update controls
load Uopt.mat
controls  = updateControls(controls, Uopt);
schedule  = updateSchedule(controls, schedule);

% FORWARD SOLVE
simResMeas = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

% save measurements information
save simResMeas.mat simResMeas;