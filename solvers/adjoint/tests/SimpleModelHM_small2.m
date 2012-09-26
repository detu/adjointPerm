% SimpleModelHM_small
% a 10 x 1 reservoir for history matching
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 10; ny = 1; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
m     = ones(nx,ny);
Kreal = ones(nx,ny);
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;

perm = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
resSolInit = initResSol(G, 0.0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('ratesMatchSmall');

% Introduce wells
radius = .1;
Wref = addWell(G, rock,[], 1         , 'Type', 'bhp' , 'Val',  25*barsa, 'Radius', radius, 'Name', 'i1');
Wref = addWell(G, rock, Wref, nx        , 'Type', 'bhp', 'Val',   20*barsa  , 'Radius', radius, 'Name', 'p1');

Wref = assembleWellSystem(G, Wref, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(Wref, 'NumSteps', 1, 'TotalTime', 0.01*day, 'Verbose', verbose);


controls = initControls(schedule, 'ControllableWells', (1:2), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% forward run
simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              
                              
 