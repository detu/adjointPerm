%compareGradients - compare gradients computed numerically and with adjoint

initSimpleModel
%{
% whether or not to show output
verbose = false;
verboseLevel = 0;

nx = 11; ny = 11; nz = 1;
G = cartGrid([nx ny nz]);
G = computeGeometry(G);

rock.perm = exp(( 3*rand(nx*ny, 1) + 1))*100;
%rock.perm(2*nx+1:2*nx+4) = 1;
rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid  = initTwoPhaseIncompressibleFluidHess('swc', 0, 'sor', 0);
resSolInit = initResSol(G, 0.0);
resSolInit.sw = .2*ones(G.cells.num, 1);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('simpleNPV');
%objectiveFunction = str2func('recovery');

% Introduce wells
radius = .1;
W = putWell(G, rock,[], 5*nx +6       , 'Type', 'rate'  ,'Val', -1  ,'Radius', radius, 'Name', 'p1');
W = putWell(G, rock, W, 1             , 'Type', 'bhp', 'Val', 10 ,'Radius', radius, 'Name', 'i1');
W = putWell(G, rock, W, nx            , 'Type', 'bhp', 'Val', 10 , 'Radius', radius, 'Name', 'i2');
W = putWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val', 10 , 'Radius', radius, 'Name', 'i3');
W = putWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val', 10  , 'Radius', radius, 'Name', 'i4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', .4*totVol, 'Verbose', verbose);

% Controlls: Only injectors are controllable and they should sum up to 1
%controls = initControls(schedule, 'ControllableWells', (2:5), ...
%                                  'RateMinMax', [0 1], ...
%                                  'LinEqConst', {ones(1,4), 1}, ...
%                                  'Verbose', verbose);
controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'RateMinMax', [0 1], ...
                                  'BHPMinMax', [100 200], ...
                                  'Verbose', verbose);

%}
% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
                              
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
grad   = computeGradient(W, adjRes, schedule, controls);  

numGrad = computeNumericalGradient(simRes, G, S, W, rock, fluid, ...
                                   schedule, controls, objectiveFunction)
adjGrad = cell2mat(grad)


figure; hold on
for k = 1 : size(numGrad, 1);
    plot(numGrad(k,:), '-ob');
    plot(adjGrad(k,:), '-xr');
end
legend('Numerical', 'Adjoint')

