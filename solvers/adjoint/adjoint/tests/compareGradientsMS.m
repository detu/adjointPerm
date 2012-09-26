% simple adjoint test
%
%clear

% whether or not to show outadd
verbose = 0;
verboseLevel = 1;

nx = 11; ny = 11; nz = 1;
G = cartGrid([nx ny nz]);
G = computeGeometry(G);

rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100*milli*darcy;
%rock.perm(2*nx+1:2*nx+4) = 1;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid  = initCoreyFluid('mu', [1 5],'sr', [0 0]);
%fluid  = initTwoPhaseIncompressibleFluid('swc', 0, 'sor', 0);
resSolInit = initResSol(G, 0.0);
resSolInit.s = .2*ones(G.cells.num, 1);

%S = assembleMimeticSystem(G, rock, 'Type', 'mixed', 'Verbose', verbose);
S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose);

% Choose objective function
%objectiveFunction = str2func('simpleNPV');
objectiveFunction = str2func('recovery');

% Introduce wells
radius = .1;
W = addWell(G, rock,[], 5*nx +6   , 'Type', 'bhp' , 'Val',  0  ,'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, 2*nx +8   , 'Type', 'bhp' , 'Val',  0  ,'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, 1         , 'Type', 'rate', 'Val',  1/4/day ,'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, nx        , 'Type', 'rate', 'Val',  1/4/day , 'Radius', radius, 'Name', 'i2');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'rate', 'Val',  1/4/day , 'Radius', radius, 'Name', 'i3');
W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp' , 'Val',  10 *barsa  , 'Radius', radius, 'Name', 'i4');
 
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', .7*totVol*day, 'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
%controls = initControls(schedule, 'ControllableWells', (2:5), ...
%                                  'RateMinMax', [0 1], ...
%                                  'LinEqConst', {ones(1,4), 1}, ...
%                                  'Verbose', verbose);
controls = initControls(schedule, 'ControllableWells', (2:6), ...
                                  'RateMinMax', [0 1]/day, ...
                                  'Verbose', 0);


% Multiscale stuff
p  = partitionCartGrid([nx ny nz], [3 3 1]);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);

p2 = partitionCartGrid([nx ny nz], [6 6 1]);
SG = generateCoarseGrid(G, p2, 'Verbose', verbose);
resSolInit.s_c = .2*ones(SG.cells.num, 1);


CS = generateCoarseSystem (G, rock, S, CG, ones(G.cells.num,1), ...
                           'Verbose', verbose, ...
                           'Overlap', 8);
                       

W = generateCoarseWellSystem(G, S, CG, CS, ones(G.cells.num,1), rock, W, ...
                             'OverlapWell', 8, ...
                             'OverlapBlock', 8);                                                 
W = addAdjointWellFields(CG, W);

CS = generateSat2MS( G, S, W, CG, CS, SG, p,fluid );                        
% forward run
simRes = runScheduleMS(resSolInit, G, CG, SG, S, CS, p, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                              
% adjoint run
adjRes = runAdjointMS(simRes, G, CG, SG, S, CS, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', verboseLevel);

grad   = computeGradient(W, adjRes, schedule, controls);  

numGrad = computeNumericalGradientMS(simRes, G, CG, SG, S, CS, p, W, rock, fluid, schedule, controls, objectiveFunction)
adjGrad = cell2mat(grad)

figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')
