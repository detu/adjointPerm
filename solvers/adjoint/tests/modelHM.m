function [W, fluid, S, resSolInit, verbose, verboseLevel] = modelHM(G, rock)
% modelHM
% a 15 x 15 reservoir for history matching
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

global numControlInterval;

nx = 15; ny = 15; nz = 1;

% Introduce wells
radius = .1;
% injectors
I   =  8 ;
J   =  8;
R   = 100/day;
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
   W = verticalWell(W, G, rock, I(i), J(i), 1, 'Type', 'rate', ...
                    'Val', -25/day , 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$']);
end


fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
fluid      = adjointFluidFields(fluid);
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
W          = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);
% welSolInit = initWellSol(W, 300*barsa());

% totVol   = sum(G.cells.volumes.*rock.poro);
% schedule = initSchedule(W, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);
% 
% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                                   'Verbose', verbose, ...
%                                   'NumControlSteps', 1);
return;
                    
 