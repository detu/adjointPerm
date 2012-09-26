% Imitation of Norne model simulation
% real geometry but 'fake' reservoir properties
% explanation can be found at
% http://www.sintef.no/Projectweb/GeoScale/Simulators/MRST/RealFieldModelII/
% 07.10.2009 - suwartad@itk.ntnu.no

addpath ../optimsecond
% load norneESeg_new;
% 
% % whether or not to show output
% verbose = false;
% verboseLevel = 0;
% 
% G = Gnew; clear Gnew;
% % G          = computeGeometry(G(1));
% K          = logNormLayers(G.cartDims, rand(9,1), 'sigma', 2);
% K          = K(G.cells.indexMap);
% K          = 200 + (K-min(K))/(max(K)-min(K))*1800;
% rock.perm  = K;                 % permeability should be multlied by *milli*darcy
% rock.poros = 0.25*(K/200).^0.1; %permeability-porosity relationship phi~0.25*K^0.11, within range 0.25-0.32
% clear K;

load NorneData

fluid  = initTwoPhaseIncompressibleFluid('muo', 5, 'muw', 1);
resSolInit = initResSol(G, 0.0);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('simpleNPV');

% display geometry of the reservoir
% clf,
% plotCellData(G,log10(rock.perm));
% axis off, view(15,60), h=colorbar('horiz');
% cs = [200 400 700 1000 1500 2000];
% caxis(log10([min(cs) max(cs)]));
% set(h, 'XTick', log10(cs), 'XTickLabel', num2str(round(cs)'));
% zoom(2.5), title('Log_{10} of x-permeability [mD]');

% introduce wells
% Set vertical injectors, completed in the lowest 12 layers.
% nz  = G.cartDims(3);
I   = [ 9, 26,  8, 25, 35, 10];
J   = [14, 14, 35, 35, 68, 75];
R   = [ 4,  4,  4,  4,  4,  4];   % should multiply by m^3/days
% I   = [25, 35, 10];
% J   = [35, 68, 75];
% R   = [ 4,  4,  4];   % should multiply by m^3/days
% nIW = 1:numel(I); 
W   = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), nz-11:nz, 'Type', 'rate', ...
                    'Val', R(i), 'Radius', 0.1, 'Name', ['I',num2str(i)], ...
                    'Comp_i', [1,0,0]);
end

% Set vertical producers, completed in the upper 14 layers
I = [17, 12, 25, 40, 15];
J = [23, 51, 51, 70, 71];
% I = 15;
% J = 71;
nPW = (1:numel(I))+max(nIW);
% nPW = (1:numel(I));
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:14, 'Type', 'bhp', ...
                    'Val', 300, 'Radius', 0.1, ...
                    'Name', ['P',num2str(i)]);
end

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Plot grid outline and the wells
clf
subplot('position',[0.02 0.02 0.96 0.96]);
plotGrid(G);
axis tight off, zoom(1.1), view(-5,58)
plotWell(G,W,'height',10);
plotGrid(G, vertcat(W(nIW).cells));
plotGrid(G, vertcat(W(nPW).cells));

% % Schedule: one time step, [0 1 PVI]
% totVol   = sum(G.cells.volumes.*rock.poros);
% schedule = initSchedule(W, 'NumSteps', 2, 'TotalTime', .5*totVol, 'Verbose', verbose);
% rateLims = [repmat([-1 -.001], 7,1); repmat([0.001 1], 4, 1)];
% 
% 
% % Controlls: Only injectors are controllable and they should sum up to 1
% controls = initControls(schedule, 'ControllableWells', (1:11), ...
%                                   'MinMax', rateLims, ...
%                                   'LinEqConst', {[ones(1,7) zeros(1,4)], -1, [zeros(1,7) ones(1,4)], 1}, ...
%                                   'NumControlSteps', 2, ...
%                                   'Verbose', verbose);
% 
% % perform derivative checks
% usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
% Uinit = [controls.well.values]';
% Uinit = Uinit(:);
% deriv_check( Uinit, 1, usr_par);
