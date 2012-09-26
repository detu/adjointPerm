% simple adjoint test
%
clear

% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 21; ny = 21; nz = 1;
G = cartGrid([nx ny nz], [5*nx 5*ny 1*nz]);
G = computeGeometry(G);

%rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100;
%rock.perm(2*nx+1:2*nx+4) = 1;
rock.perm  = ones(G.cells.num, 1)*100*milli*darcy;
rock.poro  = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu', [1 5], 'sr', [0 0]);

resSolInit = initResSol(G, 0.0);

pv = sum(poreVolume(G, rock));
totTime = 200*day;
numSteps = 5;

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose);

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

% Introduce wells
radius = .1;
W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)      , 'Type', 'bhp' , 'Val',  300*barsa  ,'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 1             , 'Type', 'rate', 'Val',  -.25*pv/totTime ,'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, nx            , 'Type', 'rate', 'Val',  -.25*pv/totTime , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, nx*ny-nx+1    , 'Type', 'rate', 'Val',  -.25*pv/totTime , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny         , 'Type', 'rate', 'Val',  -.25*pv/totTime , 'Radius', radius, 'Name', 'p4');

%W = addWell(G, rock, W, 1             , 'Type', 'bhp', 'Val',  200*barsa  , 'Radius', radius, 'Name', 'p1');
%W = addWell(G, rock, W, nx            , 'Type', 'bhp', 'Val',  200*barsa  , 'Radius', radius, 'Name', 'p2');
%W = addWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val',  200*barsa  , 'Radius', radius, 'Name', 'p3');
%W = addWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val',  200*barsa  , 'Radius', radius, 'Name', 'p4');%

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: 5 time step, [0 .7*PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', numSteps, 'TotalTime', totTime, 'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1                        
controls = initControls(schedule, 'ControllableWells', (2:5), ...
                                  'RateMinMax', [-inf 0], ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 5);                               
                               
[simRes, schedule, controls, out] = optimizeObjective(G, S, W, rock, ...
                                        fluid, resSolInit, schedule, ...
                                        controls, objectiveFunction, ...
                                        'VerboseLevel', verboseLevel);
return;                              
                              
 