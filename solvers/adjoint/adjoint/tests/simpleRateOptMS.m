% simple adjoint test
%
clear

% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 41; ny = 41; nz = 1;
G = cartGrid([nx ny nz], [5*nx 5*ny nz]);
G = computeGeometry(G);
rock.perm = ones(G.cells.num,1)*100*milli*darcy;
%rock.perm = exp(( 4*rand(nx*ny, 1) + 1))*milli*darcy;
%rock.perm(2*nx+1:2*nx+4) = 1;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu', [1 5], 'sr', [0 0]);
fluid  = adjointFluidFields(fluid);

resSolInit = initResSol(G, 0.0);

%resSolInit.s = .2*ones(G.cells.num, 1); 

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose);
S.type = 'mixed';

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

% Introduce wells
radius = .1;
W = addWell(G, rock,[], 15*nx +28 , 'Type', 'bhp', 'Val', 100*barsa, 'Radius', radius, 'Name', 'w1');
W = addWell(G, rock, W, 1         , 'Type', 'bhp', 'Val',  90*barsa, 'Radius', radius, 'Name', 'w2');
W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',  90*barsa, 'Radius', radius, 'Name', 'w3');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',  90*barsa, 'Radius', radius, 'Name', 'w4');
W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',  90*barsa, 'Radius', radius, 'Name', 'w5');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
wellSol = initWellSol(W, 0.0);
% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 10 , 'TotalTime', .5*totVol*day, ...
                        'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
% controls = initControls(schedule, 'ControllableWells', (2:5), ...
%                                   'RateMinMax', [-1 0], ...
%                                   'LinEqConst', {ones(1,4), -1}, ...
%                                   'Verbose', verbose);
                          
controls = initControls(schedule, 'ControllableWells', (2:5), ...
                                  'Verbose', verbose);
                          

% Multiscale stuff
p = partitionCartGrid([nx ny nz], [5 5 1]);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);

CS = generateCoarseSystem (G, rock, S, CG, ones(G.cells.num,1),...
                           'Verbose', verbose, ...
                           'Overlap', 4);
W = generateCoarseWellSystem(G, S, CG, CS, ones(G.cells.num,1), ...
                             rock, W, 'OverlapWell', 0, ...
                             'OverlapBlock', 4);
W = addAdjointWellFields(CG, W);                         
                         
% Make initial pressure solve:
[rS, wS] = solveIncompFlow(resSolInit, wellSol, G, S, fluid, ...
                           'wells', W, 'Solver', 'mixed');


[resSol, wellSol] = solveIncompFlowMS(resSolInit, wellSol, G, CG, p, ...
                                      S, CS, fluid, 'wells', W,      ...
                                      'Solver', 'mixed');

partNUC = partitionNonuniform(5, 25, G, rock, resSol, ...
                              'Wells',W,'wSol',wellSol, 'plotif', true);

SG =generateCoarseGrid(G,partNUC);                                 
                                
%CS = generateSat2MS( G, S, W, CG, CS, SG, p, fluid);       
resSolInit.s_c = resSolInit.s(1) * ones(SG.cells.num, 1);

[simRes1, schedule, controls, grad] = optimizeObjectiveMS(G, CG, SG, S, CS, ...
                                       p, W, rock, fluid, resSolInit, ...
                                       schedule, controls, objectiveFunction);

figure                                             
[simRes2, schedule, controls, grad] = optimizeObjective(G, S, W, rock, fluid, ...
                                                        resSolInit, schedule, ...
                                                        controls, objectiveFunction);                                             

 
