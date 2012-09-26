function normS = twophaseGravColumn(doPlot)
% Segregation flow of transport solver with gravity.
%
% Initial water saturation in all cells:
%           ______
%          |      |
%          | s=0.5|
%          |______|
%

% $Id: twophaseGravColumn.m 1829 2009-03-23 14:01:07Z bska $

if nargin == 0, doPlot = true; end


verbose = true;
if verbose,
   fprintf('**Running gravity column test of transport solvers**\n\n')
   tic
end


G = computeGeometry(cartGrid([10,1,10]), 'Verbose', verbose);

rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]);

fluid = initSimpleFluid();

% Initialize reservoir and set water saturation in all cells.
xe = initResSol(G, 0.0,[0.5,0.5])

S = assembleMimeticSystem(G, rock, 'Verbose', verbose, ...
                          'Gravity', true);
xe = solveIncompFlow(xe, [], G, S, fluid);

% WHAT??
% Hack for to avoid wrong sign on fluxes
%solution.cellFlux(:) = 0;
%solution.faceFlux(:) = 0;

%[gravity, flux, sources, porvol, dt] = ...
%   inflow_inc(G, rock, fluid, S, [], solution, []);
if doPlot
    figure;
    subplot(2,1,1)
    h = plotCellData(G, xe.s(:,1));
    view([0, 0]);
    axis tight; axis equal;

    subplot(2,1,2)
    h2 = plotCellData(G, xe.s(:,1));
    view([0, 0]);
    axis tight; axis equal;
end
% Define timesteps
tf = sum(poreVolume(G, rock)) * 100;
dt = tf / 10;

xi=xe;
%flux_impl = flux;
%verbose = false;

t = 0;
while t < tf,
   xe = explicitTransport(xe, G, [], [], dt,rock, fluid);
%{
   solution = twophaseUpwFEGrav(solution, G, timestep, sources, ...
                                flux, gravity,               ...
                                porvol, fluid,             ...
                                'verbose', verbose, 'dt', dt);
 %}
                                    
   xi = implicitTransport(xi, G, [], [], dt, rock, fluid);
%{
   solution_impl = twophaseUpwBEGrav(solution_impl, G, timestep, sources, ...
                                     flux_impl, gravity, porvol, fluid, ...
                                     'verbose', verbose);      
%}
   if max(max(abs([xe.s(:,1), xi.s(:,1)]))) > 1.00001,
      disp('ERROR: ******* Saturation >1 **********')
      break
   end
   if doPlot
      % plot saturation
      delete(h);
      subplot(2,1,1)
      h = plotCellData(G, xe.s(:,1));
      title('Water saturation - explicit solver')
      view([0, 0]); axis tight; axis equal;
      caxis([0 1]); colorbar;

      delete(h2);
      subplot(2,1,2)
      h2 = plotCellData(G, xi.s(:,1));
      title('Water saturation - implicit solver')
      view([0, 0]); axis tight; axis equal;
      caxis([0 1]); colorbar;
   end
  
   drawnow;
   t = t + dt;
end
normS = norm(xe.s(:,1) - xi.s(:,1));
