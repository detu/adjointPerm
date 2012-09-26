%
%clear
%%{
layer = 75;
satGridLims = [10 50];
CGDims      = [3 11 1];
overlap     = 10;
% whether or not to show output
verbose = 0;
verboseLevel = 0;

nx = 60; ny = 220; nz = 1;
cellDims = [nx ny nz];
physDims = [nx*20 ny*10 nz*2]*0.3048;
G = cartGrid(cellDims, physDims);
G = computeGeometry(G);

%load([SAMSIMROOTDIR, 'tests', filesep, 'data', filesep, 'Kx.mat']) 
% load data\Kx.mat
load Kx.mat

p = Kx(:,:,layer);

rock.perm = p(:)*milli*darcy;
%rock.perm = 1e-12*ones(size(p(:)));

%rock.perm = exp(( 0*rand(nx*ny, 1) + 1))*100;
%rock.perm(2*nx+1:2*nx+4) = 1;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu', [1 5], 'sr', [0 0]);
resSolInit = initResSol(G, 0.0);
resSolInit.s = .2*ones(G.cells.num, 1);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

% Introduce wells
radius = .2;
W = addWell(G, rock,[], 30 + 60*110, 'Type', 'bhp' , 'Val',  100*barsa, 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 1          , 'Type', 'bhp', 'Val',  90*barsa , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, nx         , 'Type', 'bhp', 'Val',  90*barsa , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, nx*ny-nx+1 , 'Type', 'bhp', 'Val',  90*barsa , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny      , 'Type', 'bhp', 'Val',  90*barsa , 'Radius', radius, 'Name', 'p4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 10, 'TotalTime', 1000*day, ...
                        'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
%controls = initControls(schedule, 'ControllableWells', (2:5), ...
%                                  'RateMinMax', [0 1], ...
%                                  'LinEqConst', {ones(1,4), 1}, ...
%                                  'Verbose', verbose);
controls = initControls(schedule, 'ControllableWells', (2:5), ...
                                  'RateMinMax', [-1 -.001]./day, ...
                                  'LinEqConst', {ones(1,4), -1}, ...
                                  'Verbose', verbose);

%Make initial solve:
[resSol, wellSol] = solveIncompFlow(resSolInit, [], G, S, fluid, ...
                                    'wells', W, 'Solver', 'mixed');

% Genereate non uniform coarse grid                                 
partNUC = partitionNonuniform(satGridLims(1), satGridLims(2),G,rock,resSol,...
                              'Wells',W,'wSol',wellSol, 'plotif', true);

SG =generateCoarseGrid(G,partNUC);

% compute initial function value 
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                     'VerboseLevel', verboseLevel);
        
obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
initVal   = obj.val

% Multiscale stuff
p = partitionCartGrid([nx ny nz], CGDims);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);

%p = partitionCartGrid([nx ny nz], [6 6 1]);
%SG = generateCoarseGrid(G, p, 'Verbose', verbose);

mob = ones(G.cells.num,1);

CS = generateCoarseSystem (G, rock, S, CG, mob,...
                           'Verbose', verbose, ...
                           'Overlap', overlap);
                        
                        
W = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, ...
                             'OverlapWell', overlap, ...
                             'OverlapBlock', overlap);
                          
W = addAdjointWellFields(CG, W);
                         

%CS = generateSat2MS( G, S, W, CG, CS, SG, p, fluid );       
resSolInit.s_c = resSolInit.s(1) * ones(SG.cells.num, 1);
%}                        

% %forward run
% [simRes1, schedule, controls, output1] = ...
%    optimizeObjectiveMS(G, CG, SG, S, CS, p, W, rock, fluid, resSolInit, ...
%                        schedule, controls, objectiveFunction);
figure                                             
[simRes2, schedule, controls, output2] = ...
   optimizeObjective(G, S, W, rock, fluid, resSolInit, ...
                     schedule, controls, objectiveFunction);  

                                             


