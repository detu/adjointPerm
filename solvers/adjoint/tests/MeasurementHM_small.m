% SimpleModelHM_small
% a 4 x 4 reservoir for history matching
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 4; ny = 4; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> m/K, where m is permeability multiplier
% let K = one

% m     = ones(nx,ny);
m     = [ 2 2 3 3; ...
          2 2 3 3; ...
          4 4 1 1; ...
          4 4 1 1 ];
Kreal = ones(nx,ny);
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
% K     = m ./ Kreal;
K     = Kreal ./ m;

save Kreal.mat Kreal;

% load permSmall;
% perm = reshape(permSmall', [1,G.cells.num]);
perm = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
resSolInit = initResSol(G, 0.0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('FwMatchSmall');

% Introduce wells
radius = .1;

Wref = addWell(G, rock,[], 1            , 'Type', 'bhp' , 'Val',  20.2*barsa, 'Radius', radius, 'Name', 'i1');
Wref = addWell(G, rock, Wref, nx        , 'Type', 'bhp', 'Val',   20*barsa  , 'Radius', radius, 'Name', 'p1');
Wref = addWell(G, rock, Wref, nx*ny-nx+1, 'Type', 'bhp', 'Val',   20*barsa  , 'Radius', radius, 'Name', 'p3');
Wref = addWell(G, rock, Wref, nx*ny     , 'Type', 'bhp', 'Val',   20*barsa  , 'Radius', radius, 'Name', 'p3');


Wref = assembleWellSystem(G, Wref, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(Wref, 'NumSteps', 4, 'TotalTime', 0.01*day, 'Verbose', verbose);


% controls = initControls(schedule, 'ControllableWells', [], ...
%                                   'Verbose', verbose, ...
%                                   'NumControlSteps', 1);
% forward run
simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              
                              
 