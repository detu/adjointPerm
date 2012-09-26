% Estimation permeability using 101 realizations 
% the 101st realization is assumed as the truth realization
% 45x45 reservoir, 10 wells: 5 injectors at one side and 5 producers at the
% other side
% Measurements: water production at the wells
% This code is to run production optimization to get optimal well rates and
% to be used in history matching, also record water cut 

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
% injector wells
W = addWell(G, rock,W, 1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', 'I1');
for i=1:5
    W = addWell(G, rock,W, (9*i - 1)*nx +1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', ['I',num2str(i+1)]);
end
% producer wells
W = addWell(G, rock,W, 45     , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', 'P1');
for i=1:5
    W = addWell(G, rock,W, i*9*nx      , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', ['P',num2str(i+1)]);
end

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol     = sum(G.cells.volumes.*rock.poro);
schedule   = initSchedule(W, 'TimeSteps', 120*day:120*day:2400*day, 'Verbose', verbose);

injMaxMin  = repmat( [1e-15 1e3]*meter^3/day, 6, 1);
prodMaxMin = repmat( [-1e3 -1e-15]*meter^3/day, 6, 1);

controls   = initControls(schedule, 'ControllableWells', (1:12), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);
                   
%---OPTIMIZATION PART---------                   
objectiveFunction = str2func('simpleNPV');
% initial  control input
U       = [controls.well.values]';
U       = U(:);
minMax  = vertcat(controls.well.minMax);
numstep = numel(schedule);
minMax  = repmat(minMax,numstep,1);
lb      = minMax(:,1);
ub      = minMax(:,2);

auxdata   = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};

% Specify some model characteristics.
A   = [];
b   = [];

% construct Aeq, equality matrix constraint
Aeq = zeros(20,120);
for k=1:20
    startInd = 12*k-11;
    endInd   = 12*k;
    Aeq(k,startInd:endInd) = ones(1,12);
end

beq = zeros(20,1);

% run KNITRO
tic;
[U,fval,exitflag,output,lambda] = ktrlink(@(U)NPVfunction(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],[], 'knitro.opt');
toc;

% save optimal well rates
Uopt = U;
Uopt = convertTo(Uopt,meter^3/day);
save Uopt.mat Uopt;
