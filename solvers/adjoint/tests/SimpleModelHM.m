% SimpleModelHM
% a 15 x 15 reservoir for history matching
%
% whether or not to show output
clear; clc;
verbose      = false;
verboseLevel = 0;

nx = 15; ny = 15; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

load perm;

m     = reshape(m, [1,G.cells.num]);
Kreal = rot90(perm);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
save Kmodel.mat K;

perm = reshape(K, [1,G.cells.num]);

perm = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);


% Choose objective function
% objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');
% objectiveFunction = str2func('Jac_Swn');

% Introduce wells
radius = .1;
% injectors
% I   = [ 1, 1, 1 , 1 , 6, 6, 6 , 6 , 10, 10, 10, 10, 15, 15, 15, 15];
% J   = [ 1, 5, 11, 15, 1, 5, 11, 15, 1 , 5 , 11, 15, 1 ,  5, 11, 15];
% R   = 4*[ 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1 , 1, 1, 1 ,  1, 1, 1]*1000*meter^3/day;
I   =  8 ;
J   =  8;
R   = 4*1000*meter^3/day;
nIW = 1:numel(I); W = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), nz:nz, 'Type', 'rate', ...
                    'Val', R(i), 'Radius', 0.1, 'Comp_i', [1,0,0], ...
                    'name', ['I$_{', int2str(i), '}$']);
end

% producers
I   = [ 3, 3 , 12, 12 ];
J   = [ 4, 13,  4 , 13 ];
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1, 'Type', 'bhp', ...
                    'Val', 350*barsa(), 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$']);
end


fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
W          = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 0.0);
% welSolInit = initWellSol(W, 300*barsa());

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 360*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% forward run
simRes_ref = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                          
return;                              
                              
 