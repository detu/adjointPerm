%function solution = twophaseQuarterFiveSpotExp(dims)
% simple example of cartesian grid with unit permeability
% compare implicit and explicit solver of saturation equation

dims=[10 10 10];

verbose = true;

% To export this function's workspace to console, uncomment line below
%c = onCleanup(@()exportWorkspace);

G = cartGrid(dims);
G = computeGeometry(G, 'Verbose', verbose);

rock.poro = repmat(0.3,[G.cells.num, 1]);
rock.perm  = repmat(100,[G.cells.num, 1]);

fluid = initSimpleFluid();

solution = initResSol(G, 0.0);

S = assembleMimeticSystem(G, rock, 'Verbose', verbose, 'Gravity', false);

someCells = [1;G.cells.num];
sources   = [1;-1];
src       = addSource([], someCells, sources, [1,0,0;1,0,0]);

solution     = solveIncompFlow(solution, [], G, S, fluid, 'src', src);
solution_exp = solution;
[gm, sources, porvol] = inflow(G, rock, S, [], solution, [], 'src', src);


%% plot output
figure;
clf
subplot(2,2,1)
   plotGrid(G, [], 'FaceColor', 'none');
   plotGrid(G, someCells, 'FaceColor', 'r');
   title('Source-cells')
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;

subplot(2,2,2)
   plotCellData(G, solution.cellPressure);
   title('Pressure')
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;

subplot(2,2,3)
   plotCellData(G, sqrt(S.C'*abs(solution.cellFlux)) );
   title('Sqrt - Flux intensity')
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;

subplot(2,2,4)
   h = plotCellData(G, solution.s(:,1));
   title('Initial water saturation')
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;

figure;
subplot(1,2,1)
   h = plotCellData(G, solution.s(:,1));
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;
subplot(1,2,2)
   h2 = plotCellData(G, solution_exp.s(:,1));
   view(3); camproj perspective; axis tight; axis equal; camlight headlight;

% define timesteps
tf = sum(porvol) / 2;
timestep = tf / 10;
dt = min(timestep, 0.1);
t  = 0;

while t < tf,
   solution     = twophaseUpwBE(solution, timestep, sources, gm, ...
                                porvol, fluid, 'verbose', verbose, ...
                                'nltol', 1e-8);

   solution_exp = twophaseUpwFE(solution_exp, timestep, sources, gm, ...
                                porvol, fluid, 'verbose', verbose,   ...
                                'dt', dt);

   subplot(2,1,1)
      h = plotCellData(G, solution.s(:,1));
      title('Water saturation - implicit solver')
      view(3);
      axis tight; axis equal;
      caxis([0 1]); colorbar;

   subplot(2,1,2)
      h2 = plotCellData(G, solution_exp.s(:,1));
      title('Water saturation - explicit solver')
      view(3);
      axis tight; axis equal;
      caxis([0 1]); colorbar;

   drawnow;
   t = t + dt;
end
