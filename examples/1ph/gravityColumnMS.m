%% Multiscale Pressure Solver: Simple Case Driven by Gravity
% Compare the fine-scale and the multiscale solver for the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad
%    v=\textbf{--}\frac{K}{\mu} \bigl[\nabla p+\rho g\nabla z\bigr],$$
%
% This example is a direction continuation of <gravityColumn.html "My First
% Flow-Solver: Gravity Column"> and introduces the multiscale flow solver
% without going into specific details. More details can be found in the
% <simpleBCMS.html "Basic Multiscale Tutorial">.

%% Define the model
% The domain [0,1]x[0,1]x[0,30] is discretized using a Cartesian grid with
% homogeneous isotropic permeability of 100 mD. The fluid has density 1000
% kg/m^3 and viscosity 1 cP and the pressure is 100 bar at the top of the
% structure
nx = 2; ny = 2; nz = 30;
Nx = 1; Ny = 1; Nz =  6;
gravity reset on
G          = cartGrid([nx, ny, nz], [1, 1, 30]);
G          = computeGeometry(G);
rock.perm  = repmat(0.1*darcy(), [G.cells.num, 1]);
fluid      = initSingleFluid();
bc         = pside([], G, 'TOP', 1:nx, 1:ny, 100.*barsa());

%% Assemble and solve the fine-scale linear system
S   = computeMimeticIP(G, rock);
sol = solveIncompFlow(initResSol(G , 0.0), initWellSol([], 0.0), ...
                      G, S, fluid,'bc', bc);

%% Plot the fine-scale solution
newplot;
subplot(3, 2, [1 3])
   plotFaces(G, 1:G.faces.num, sol.facePressure ./barsa());
   set(gca, 'ZDir', 'reverse'), title('Fine-scale pressure [bar]')
   view(45,5), cx = caxis; colorbar

%% Multiscale system
p  = partitionUI(G, [Nx, Ny, Nz]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);

CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]),'bc', bc);
xrMs = solveIncompFlowMS (initResSol(G, 0.0), [], G, CG, p, ...
                          S, CS, fluid, 'bc', bc);

%% Plot the coarse-scale solution
% As one clearly can see from the plot, the multiscale solution only
% captures the gravity effect on the coarse scale. To also capture
% fine-scale gravity effects, one can either add extra correction functions
% or insert the multiscale solution into the fine-scale equations and solve
% for a residual correction
subplot(3, 2, [2 4])
   plotFaces(G, 1:G.faces.num, xrMs.facePressure./barsa());
   set(gca, 'ZDir', 'reverse'); title('Coarse-scale pressure [bar]')
   view(45,5); caxis(cx); colorbar
subplot(3, 2, [5 6]);
   plot(1:nz, sol.cellPressure(1:nx*ny:nx*ny*nz) ./barsa(),'-o',...
        1:nz, xrMs.cellPressure(1:nx*ny:nx*ny*nz)./barsa(),'-*');
   legend('fine','coarse',2);

%%
% #COPYRIGHT_EXAMPLE#

%%
% <html>
% <font size="-1">
%   Last time modified: 
%   $Id: gravityColumn.m 1907 2009-03-26 14:40:38Z knl $
% </font>
% </html>
displayEndOfDemoMessage(mfilename)
