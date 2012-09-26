function twophaseGravVertical(doPlot)

% simple example of cartesian grid with unit permeability
% using an explicit transport solver with gravitation.
% Initial water saturation is:
%       _________________
%      |        |        |
%      | s1 = 1 | s2 = 0 |
%      |________|________|

% $Id: twophaseGravVertical.m 1829 2009-03-23 14:01:07Z bska $

if nargin == 0
    doPlot = true;
end

dims = [10, 1, 10];

verbose = true;
if verbose,
   fprintf(['** Running test of transport solvers with gravity - ', ...
            'initial vertical split saturation **\n\n'])
end

G = cartGrid(dims);
G = computeGeometry(G, 'Verbose', verbose);
% G.faces.normals = G.faces.normals * -1;
rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]);

fluid = initSimpleFluid();
resSol = initResSol(G, 0.0);
resSol.s(:,1) = repmat([0.8*ones(dims(1)/2,1);  ...
                        0.0*ones(dims(1)/2,1)], [dims(3), 1]);

S = assembleMimeticSystem(G, rock, 'Verbose', false, 'Gravity', true);
resSol = solveIncompFlow(resSol, [], G, S, fluid);

if doPlot,
   % Find internal faces for plotting facefluxes
   intFInx = find(all(G.faces.neighbors > 0, 2));

   figure;
   subplot(2,2,1),
      plotCellData(G, resSol.cellPressure);
      title('Pressure')
      view([0, 0]);
      axis tight; axis equal;

   subplot(2,2,2),
      plotCellData(G, sqrt(S.C'*abs(resSol.cellFlux)) );
      title('Sqrt - Flux intensity')
      view([0, 0]);
      axis tight; axis equal;

   subplot(2,2,3),
      newplot();
      plotFaces(G, intFInx, resSol.faceFlux(intFInx));
      colorbar; title('Faceflux')
      view(3);
      axis tight; axis equal;

   subplot(2,2,4),
      plotCellData(G, resSol.s(:,1));
      title('Initial water saturation')
      view([0, 0]);
      axis tight; axis equal;

   satf = figure();
   subplot(2,1,1),
      h = plotCellData(G, resSol.s(:,1));
      view([0, 0]);
      axis tight; axis equal; colorbar;

   subplot(2,1,2),
      h2 = plotCellData(G, resSol.s(:,1));
      view([0, 0]);
      axis tight; axis equal;
      caxis([0 1]);
end
% Define timesteps
tf       = 3000;
timestep = 100;
t        = 0;

resSol_impl = resSol;
verbose = false;

while t < tf,
   resSol = explicitTransport(resSol, G, [], [], timestep, ...
                              rock, fluid, 'verbose', verbose);

   resSol_impl = implicitTransport(resSol_impl, G, [], [], timestep, ...
                                   rock, fluid, 'verbose', verbose);

   if max(max(abs([resSol.s(:,1), resSol_impl.s(:,1)]))) > 1.00001,
      disp('ERROR: ******* Saturation exceeds 1 **********')
      break
   end

   if doPlot,
      % plot saturation
      figure(satf)
      delete(h);
      subplot(2,1,1),
         h = plotCellData(G, resSol.s(:,1));
         axis tight, axis equal
         title('Water saturation - explicit solver')
         view([0, 0]);
         caxis([0, 1]);
         colorbar;

      figure(satf)
      delete(h2);
      subplot(2,1,2),
         h2 = plotCellData(G, resSol_impl.s(:,1));
         title('Water saturation - implicit solver')
         view([0, 0]);
         axis tight, axis equal
         caxis([0, 1]);
         colorbar
   end

   % Update solution of pressure equation.
   resSol = solveIncompFlow(resSol, [], G, S, fluid);
   resSol_impl = solveIncompFlow(resSol_impl, [], G, S, fluid);

   drawnow;
   t = t + timestep;
end

normSol = num2str(norm(resSol.s(:,1) - resSol_impl.s(:,1)));
