function normS = twophaseQuarterFiveSpotGrav(doPlot)
% Simple example of cartesian grid with unit permeability
% Compare explicit solvers of saturation equation with and without gravity.
% (compare twophaseUpwFEGrav against twophaseUpwFE)

if nargin == 0
    doPlot = true;
end

dims = [10, 1, 10];

verbose = true;
if verbose,
   fprintf(['** Running quarter 5-spot test of transport ',...
            'solvers with gravity **\n\n']);
   tic
end

G = cartGrid(dims);
G = computeGeometry(G, 'Verbose', verbose);
%G.faces.normals = G.faces.normals*(-1);


rock.poro = repmat(0.3,[G.cells.num, 1]);
rock.perm  = repmat(100,[G.cells.num, 1]);

fluid = initSimpleFluid();

resSol = initResSol(G, 0.0);

S = assembleMimeticSystem(G, rock, 'Verbose', verbose, 'Gravity', false);
% Test solvers are without gravity:
someCells = [1;G.cells.num];
sources   = [ 1; -1];
src       = addSource([], someCells, sources, [1,0,0; 1,0,0]);

resSol = solveIncompFlow(resSol, [], G, S, fluid, 'src', src);

% Increase gravity effects
% S.constants.gravity = S.constants.gravity * 100;

%[gm] = inflow(G, rock, S, [], resSol, [], 'src', src);
[grav, flux, sources, porvol] = ...
   inflow_inc(G, rock, fluid, S, [], resSol, [], ...
              'src', src, 'ComputeDt', true);


%% Plot output
if doPlot
   figure;
   subplot(2,1,1),
      h = plotCellData(G, resSol.s(:,1));
      view([0, 0]), axis tight equal

   subplot(2,1,2),
      h2 = plotCellData(G, resSol.s(:,1));
      view([0, 0]), axis tight equal
end
% Define timesteps
tf = sum(porvol)/2;
timestep = sum(porvol)/10;
t  = 0;

resSol_impl = resSol;
flux_impl = flux;
verbose = true;
dt = 0.0504;
while t < tf,
   resSol      = twophaseUpwFEGrav(resSol, G, timestep, sources, ...
                                   flux, grav, porvol, fluid, ...
                                   'verbose', verbose, 'dt', dt);

   resSol_impl = twophaseUpwBEGrav(resSol_impl, G, timestep, sources, ...
                                   flux_impl, grav, porvol, fluid,    ...
                                   'verbose', verbose);

   if max(max(abs(resSol.s(:,1))), max(abs(resSol_impl.s(:,1)))) > 1+1.0e-5
      disp('ERROR: ******* Saturation exceeds 1 **********')
      break
   end

   % Plot saturation
   if doPlot,
      delete(h);
      subplot(2,1,1),
         h = plotCellData(G, resSol.s(:,1));
         title('Water saturation - explicit solver')
         view([0, 0]), axis tight equal
         caxis([0 1]), colorbar

      delete(h2);
      subplot(2,1,2),
         h2 = plotCellData(G, resSol_impl.s(:,1));
         title('Water saturation - implicit solver')
         view([0, 0]), axis tight equal
         caxis([0 1]), colorbar

      drawnow
   end

   % Update solutions of pressure equation.
   resSol = solveIncompFlow(resSol, [], G, S, fluid, 'src', src);

   [grav, flux, sources, porvol] = inflow_inc(G, rock, fluid, S, ...
                                              [], resSol, [],    ...
                                              'src', src,        ...
                                              'ComputeDt', true);

   resSol_impl = solveIncompFlow(resSol_impl, [], G, S, fluid, ...
                                 'src', src);

   [grav, flux_impl] = inflow_inc(G, rock, fluid, S, [], ...
                                  resSol_impl, [],       ...
                                  'src', src, 'ComputeDt', true);

   t = t + timestep;
end

normS =   num2str(norm(resSol.s(:,1) - resSol_impl.s(:,1)));
