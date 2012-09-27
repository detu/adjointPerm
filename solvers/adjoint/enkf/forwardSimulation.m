function [X outputData] = forwardSimulation(X0, U, auxdata)


%loading parameters
[G, rock, fluid,wellPar,simParam,modelPar] = deal(auxdata{:});

verbose = modelPar.verbose;
verboseLevel = modelPar.verboseLevel;

wellRadius = wellPar.wellRadius;
prodIdx = wellPar.prodIdx;
injIdx = wellPar.injIdx;
objectiveFunction= simParam.objectiveFunction;
simTime = simParam.simTime;
numContrSteps = simParam.numContrSteps;
totRate = simParam.totRate;

% Introduce wells
% set zero-poros cells inactive
numCells = G.cells.num;
numProd = length(prodIdx);

rock.perm = X0((numProd+1):(numProd +numCells),1); % from the state %pReal(:,trueRealization)*darcy;

% To fix
S01 = X0( (numProd +numCells+1): (numProd +2*numCells) , 1);
P01 = X0( (numProd +2*numCells+1): (numProd +3*numCells) , 1);

% In this way we start always from the same initial conditions
S0 = 0.1;
P0 = 0;

resSolInit = initResSol(G, P0, S0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_simple');

radius = wellRadius;
numInj = length(injIdx);
numProd = length(prodIdx);

for j1=1:numInj
    W = addWell(G, rock,[], injIdx(j1), 'Type', 'rate' ,'Val',     totRate , 'Radius', radius, 'Name', ['i' num2str(j1)]);
end

for j1=1:numProd
    W = addWell(G, rock,W, prodIdx(j1), 'Type', 'rate' ,'Val',   -totRate/numProd , 'Radius', radius, 'Name', ['p' num2str(j1)]);
end


W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');


% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', numContrSteps, 'TotalTime', simTime*day, 'Verbose', verbose);


injMaxMin  = repmat( [0.001 1000000]*meter^3/day, 1);
%injMaxMin  = repmat( [6e2 1000000]*meter^3/day, 1);
prodMaxMin = repmat( [-1000000 -0.001]*meter^3/day, 4, 1);

% injMaxMin  = repmat( [1e3 6e3]*meter^3/day, 1);
% prodMaxMin = repmat( [-1e4 -0.001]*meter^3/day, 4, 1);

controls = initControls(schedule, 'ControllableWells', (1:(numInj+numProd)), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', {ones(1,5), 0}, ...
                       'Verbose', verbose);

                   % initial  control input
                   U = [controls.well.values]';
                   U = U(:);
                   
                   % update schedule
                   controls  = updateControls(controls, U);
                   schedule  = updateSchedule(controls, schedule);






% FORWARD SOLVE
simRes    = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

% COMPUTE NPV
obj = objectiveFunction(G, S, W, rock, fluid, simRes);

% % plot Watercut
% fwProd = waterCut(G, S, W, rock, fluid, simRes);
% figure(44)
% clf
% hold on
% plot(1:numContrSteps,ones(1,numContrSteps)*(80/(80+20)),'--k')
% legend('Economic limit')
% hold off
% hold on
% plot(fwProd')
% hold off
% 
% %% plot saturations and pressures at end time
% resSol  = simRes(end).resSol;
% %wellSats  = resSol.s( wellCells );
% Sats  = resSol.s;
% Ps =  resSol.cellPressure;
% figure(45);
% clf
% plotCellData(G,Sats); axis square; axis tight; colorbar;
% title('Water saturation')
% figure(46);
% clf
% plotCellData(G,convertTo(Ps, barsa())); axis square; axis tight; colorbar;
% title('Pressure (bar)')

f_RO   = -obj.val/1e7 %normalization

%fprintf('Current NPV: %e\n', obj.val);

if nargout > 1
    outputData = {simRes,G};   
end



X = [fwProd(:,end); rock.perm; resSol.s; resSol.cellPressure];