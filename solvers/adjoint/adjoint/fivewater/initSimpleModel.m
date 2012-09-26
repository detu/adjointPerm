% InitSimpleModel
%
% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 10; ny = 10; nz = 1;
G = cartGrid([nx ny nz]);
G = computeGeometry(G);

rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100;
% rock.perm = ones( size(rock.perm) );
rock.poros = repmat(0.3, [G.cells.num, 1]);

fluid  = initTwoPhaseIncompressibleFluid2('muo', 5, 'muw', 1);
resSolInit = initResSol(G, 0.0);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose);

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

% Introduce wells
radius = .1;
W = putWell(G, rock,[], 1           , 'Type', 'bhp', 'Val',   6  , 'Radius', radius, 'Name', 'i1');
W = putWell(G, rock, W, nx          , 'Type', 'bhp', 'Val',   5  , 'Radius', radius, 'Name', 'p1');
W = putWell(G, rock, W, nx*ny-nx+1  , 'Type', 'bhp', 'Val',   5   , 'Radius', radius, 'Name', 'p3');
W = putWell(G, rock, W, nx*ny       , 'Type', 'bhp', 'Val',   5   , 'Radius', radius, 'Name', 'p3');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poros);
schedule = initSchedule(W, 'NumSteps', 10, 'TotalTime', totVol, 'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
controls = initControls(schedule, 'ControllableWells', (1:4), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 10);
                          
return;                              
                              
 