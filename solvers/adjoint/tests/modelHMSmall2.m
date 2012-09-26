% SimpleModelHM_small
% a 10 x 1 reservoir for history matching
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
% fluid      = initSingleFluid('mu', [1 5], 'sr', [0 0]);
resSolInit = initResSol(G, 0.0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('ratesMatchSmall');

% Introduce wells
radius = .1;
W = addWell(G, rock,[], 1         , 'Type', 'bhp' , 'Val',  25*barsa, 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',   20*barsa  , 'Radius', radius, 'Name', 'p1');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', .01*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:2), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);