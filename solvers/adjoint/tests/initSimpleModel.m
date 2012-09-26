% InitSimpleModel
%
% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 4; ny = 4; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100*milli*darcy;
rock.perm = ones( size(rock.perm) )*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
resSolInit = initResSol(G, 0.0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');
% objectiveFunction = str2func('Jac_Swn');

% Introduce wells
radius = .1;
% W = addWell(G, rock,[], 1         , 'Type', 'rate' , 'Val', 1.5/day, 'Radius', radius, 'Name', 'i1');
% W = addWell(G, rock, W, nx        , 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p1');
% W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p3');
% W = addWell(G, rock, W, nx*ny     , 'Type', 'rate', 'Val',  -.5/day  , 'Radius', radius, 'Name', 'p3');
% W = addWell(G, rock,[], 1         , 'Type', 'bhp' , 'Val',  1e-5*barsa, 'Radius', radius, 'Name', 'i1');
% W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',   .2e-5*barsa  , 'Radius', radius, 'Name', 'p1');
% W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',   .2e-5*barsa  , 'Radius', radius, 'Name', 'p3');
% W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',   .2e-5*barsa  , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock,[], 1         , 'Type', 'bhp' , 'Val',  1, 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',   .2  , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',   .2  , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',   .2  , 'Radius', radius, 'Name', 'p3');


W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 2, 'TotalTime', totVol*day, 'Verbose', verbose);


controls = initControls(schedule, 'ControllableWells', (1:4), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 2);
                          
return;                              
                              
 