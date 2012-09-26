% silvias case
%
clear

% whether or not to show output
verbose = false;
verboseLevel = 0;

% load([SAMSIMROOTDIR, 'tests', filesep, 'data', filesep, 'silvia.mat']) 
load silvia.mat; 

nx = 31; ny = 41; nz = 1;
cellDims = [nx ny nz];
physDims = [100*nx 100*ny 100*nz]*0.3048;

G = cartGrid(cellDims, physDims);

% set zero-poros cells inactive
rock.poro = max(rock.poro, 1e-3);

G = computeGeometry(G);

fluid  = initCoreyFluid('mu', [0.3 3], 'sr', [.2 .2]);
resSolInit   = initResSol(G, 0.0, 0.2);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_simple');

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

totRate = 5000*0.1590/day; % cubic meters
% Introduce wells
radius = .1;
W = addWell(G, rock,[], nx*20 + 16, 'Type', 'rate' ,'Val',     totRate , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, 1         , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, nx        , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i2');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i3');
W = addWell(G, rock, W, nx*ny     , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', 1000*day, 'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
controls = initControls(schedule, 'ControllableWells', (2:5), ...
                                  'RateMinMax', [-totRate -0.001/day], ...
                                  'LinEqConst', {ones(1,4), -totRate}, ...
                                  'Verbose', verbose);

% injMaxMin  = repmat( [0.001 1000000]*meter^3/day, 1);
% prodMaxMin = repmat( [-1000000 -0.001]*meter^3/day, 4, 1);
% 
% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                        'MinMax', [injMaxMin; prodMaxMin],...
%                        'LinEqConst', {ones(1,5), 0}, ...
%                        'Verbose', verbose);
                          
% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                                   'RateMinMax', [-totRate -0.001/day], ...
%                                   'LinEqConst', {ones(1,5), 0}, ...
%                                   'Verbose', verbose);
                            
[simRes, schedule, controls, grad] = optimizeObjective(G, S, W, rock,     ...
                                                       fluid, resSolInit, ...
                                                       schedule, controls,...
                                                       objectiveFunction);
return;                              
                              
 