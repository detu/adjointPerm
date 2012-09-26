%% Two-Point Flux Approximation Solvers with gravity
gravity reset on


%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured format in which cells, faces, nodes, etc. are given
% explicitly.
nx = 20; ny = 20; nz = 20;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G, 'Verbose', true);

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
rock.perm = repmat(100, [G.cells.num, 1]) .* 1e-3*darcy();
fluid     = initSingleFluid();

%% Introduce wells
% We will include two vertical pressure-controlled wells. The wells are described
% using a Peacemann model, giving an extra set of (trivial) equations that need to be
% assembled. We need to specify ('InnerProduct', 'ip_tpf') to get the
% correct well model for TPFA.
%
cellsWell1 =  1 : nx*ny : nx*ny*nz;
W = addWell(G, rock, [], cellsWell1,          ...
            'Type', 'bhp', 'Val', 2.2049*barsa(), ...
            'InnerProduct', 'ip_tpf');

cellsWell2 = nx*ny: nx*ny : nx*ny*nz;
W = addWell(G, rock,  W, cellsWell2,      ...
            'Type', 'bhp' , 'Val', 1.0*barsa(), ...
            'InnerProduct', 'ip_tpf');


%% APPROACH 1: Direct/Classic TPFA
% Initialize solution structures for reservoir and wells.
resSol1  = initResSol(G, 0.0);
wellSol1 = initWellSol(W, 1.0*barsa());

% Compute one-sided transmissibilities.
T = computeTrans(G, rock, 'Verbose', true);

% Solve linear system construced from T and W to obtain solution for flow
% and pressure in the reservoir and the wells. Notice that the TPFA solver
% is different from the one used for mimetic systems.
[resSol1, wellSol1] = incompTPFA(resSol1, wellSol1,  G, T, ...
                                 fluid, 'wells', W, 'Verbose', true);


%% APPROACH 2: Mimetic with TPFA-inner product
% Initialize solution structures for reservoir and wells.
resSol2 = initResSol(G, 0.0);
wellSol2 = initWellSol(W, 1.0*barsa());

% Compute mimetic innerproduct equivalent to two-point flux for Cartesian
% grids.
IP = computeMimeticIP(G, rock, 'Verbose', true, ...
                      'InnerProduct', 'ip_tpf');

%%
% Solve mimetic linear hybrid system
[resSol2, wellSol2] = solveIncompFlow(resSol2, wellSol2, ...
                                     G, IP, fluid, 'wells', W);

%% Report results
% Report pressure drop computed by the two solvers.
dP1 = (wellSol1(1).pressure-wellSol1(2).pressure) ./ barsa();
dP2 = (wellSol2(1).pressure-wellSol2(2).pressure) ./ barsa();
disp(['DeltaP,  direct TPFA: ', num2str(dP1)])
disp(['DeltaP, mimetic TPFA: ', num2str(dP2)])

%%
% Plot the pressure and producer inflow profile
clf
subplot('Position', [0.05,0.55,0.4, 0.35])
   plotCellData(G, resSol1.cellPressure ./ barsa());
   title('Pressure: direct TPFA with pressure control')
   view(45, 25), camproj perspective, axis tight off, camlight headlight
   cax = caxis;

subplot('Position', [0.55,0.55,0.4, 0.35])

   plotCellData(G, resSol2.cellPressure ./ barsa());
   title('Pressure: mimetic TPFA with pressure control')
   view(45, 25), camproj perspective, axis tight off, camlight headlight
   caxis(cax);
   
subplot('Position', [0.15,0.4,0.7, 0.08])

colorbar south; caxis(cax);axis tight off;
   
subplot('position', [0.1, 0.1, 0.8, 0.25])
   plot(-wellSol1(2).flux .* day(), 'b-*'); hold on
   plot(-wellSol2(2).flux .* day(), 'r--');
   legend('Direct','Mimetic')
   title('Producer inflow profile [m^3/d]');


%% Rate controlled wells
W = addWell(G, rock, [], cellsWell1,          ...
            'Type', 'RaTe', 'Val', 5.0/day(), ...
            'InnerProduct', 'ip_tpf');
W = addWell(G, rock,  W, cellsWell2,      ...
            'Type', 'rAtE' , 'Val', -5.0/day(), ...
            'InnerProduct', 'ip_tpf');


%% APPROACH 1: Direct/Classic TPFA
resSol1  = initResSol (G, 0.0);
wellSol1 = initWellSol(W, 1.0*barsa());

[resSol1, wellSol1] = incompTPFA(resSol1, wellSol1,  G, T, ...
                                 fluid, 'wells', W, 'Verbose', true);


%% APPROACH 2: Mimetic with TPFA-inner product
resSol2  = initResSol (G, 0.0);
wellSol2 = initWellSol(W, 1.0*barsa());

[resSol2, wellSol2] = solveIncompFlow(resSol2, wellSol2, ...
                                     G, IP, fluid, 'wells', W);

%% Report results
% Report pressure drop computed by the two solvers.
dP1 = (wellSol1(1).pressure-wellSol1(2).pressure) ./ barsa();
dP2 = (wellSol2(1).pressure-wellSol2(2).pressure) ./ barsa();
disp(['DeltaP,  direct TPFA: ', num2str(dP1)])
disp(['DeltaP, mimetic TPFA: ', num2str(dP2)])

%%
% Plot the pressure and producer inflow profile
clf
subplot('Position', [0.05,0.55,0.4, 0.35])
   plotCellData(G, resSol1.cellPressure ./ barsa());
   title('Pressure: direct TPFA with rate control')
   view(45, 25), camproj perspective, axis tight off, camlight headlight
   cax = caxis;

subplot('Position', [0.55,0.55,0.4, 0.35])

   plotCellData(G, resSol2.cellPressure ./ barsa());
   title('Pressure: mimetic TPFA with rate control')
   view(45, 25), camproj perspective, axis tight off, camlight headlight
   caxis(cax);
   
subplot('Position', [0.15,0.4,0.7, 0.08])

colorbar south; caxis(cax);axis tight off;
   
subplot('position', [0.1, 0.1, 0.8, 0.25])
   plot(-wellSol1(2).flux .* day(), 'b-*'); hold on
   plot(-wellSol2(2).flux .* day(), 'r--');
   legend('Direct','Mimetic')
   title('Producer inflow profile [m^3/d]');

   
%%
% #COPYRIGHT_EXAMPLE#

%%
% <html>
% <font size="-1">
%   Last time modified:
%   $Id: secondTPFA.m 2058 2009-04-21 13:25:31Z jrn $
% </font>
% </html>
displayEndOfDemoMessage(mfilename)

