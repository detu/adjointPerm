% EnKF implementation of 45x45 reservoir
% 12 wells: 6 injectors and 6 producers
% Measurement data are water-cut at each well


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
load permRealizations.mat;

% reservoir setting
nx = 21; ny = 21; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

% rock properties
m     = reshape(m, [1,G.cells.num]);
Kreal = perm(:,76);    %truth permeability no. 76 from realization
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

% change permeability realizations
permeability       = perm;
permeability(:,76) = [];

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
W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'rate' ,'Val',   totRate  , 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 5*nx+6        , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, 6*nx-5        , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, 16*nx+5       , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, 17*nx-5       , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol    = sum(G.cells.volumes.*rock.poro);
schedule  = initSchedule(W, 'NumSteps', 20, 'TotalTime', 100*totVol, 'Verbose', verbose);
%% ------------OPTIMIZATION PART-------------------------------------
injMaxMin  = repmat( [1e-15 2.5]*meter^3/day, 1, 1);
prodMaxMin = repmat( [-1 -1e-15]*meter^3/day, 4, 1);

controls   = initControls(schedule, 'ControllableWells', (1:5), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);

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
Aeq = zeros(20,100);
for k=1:20
    startInd = 5*k-4;
    endInd   = 5*k;
    Aeq(k,startInd:endInd) = ones(1,5);
end

beq = zeros(20,1);

% run KNITRO
tic;
[U,fval,exitflag,output,lambda] = ktrlink(@(U)NPVfunction(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],[], 'knitro.opt');
toc;

% save optimal well rates
Uopt = U;
Uopt = convertTo(Uopt,meter^3/day);
save UoptCase4.mat Uopt;
