function solution = twophaseQuarterFiveSpot(varargin)
% simple example of cartesian grid with unit permeability
if nargin == 1
  dims = varargin{:};
else
  dims = [30,30,1];
end

verbose = false;

% To export this function's workspace to console, uncomment line below
c = onCleanup(@()exportWorkspace);

G = cartGrid(dims);
G = computeGeometry(G, 'Verbose', verbose);

rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]);

fluid = initSimpleFluid();

solution = initResSol(G, 0.0);

S = assembleMimeticSystem(G, rock, 'Verbose', verbose, 'Gravity', false);

someCells = [1; G.cells.num];
sources   = [1; -1];
src       = addSource([], someCells, sources, [1,0,0; 1,0,0]);
solution  = solveIncompFlow(solution, [], G, S, fluid, 'src', src);

[gm, sources, porvol] = inflow(G, rock, S, [], solution, [], 'src', src);


%% plot output
clf
subplot(2,2,1)
   plotGrid(G, 'FaceColor', 'none');
   plotGrid(G, someCells, 'FaceColor', 'r');
   title('Source-cells')
   set_view()

subplot(2,2,2)
   plotCellData(G, solution.cellPressure);
   title('Pressure')
   set_view()

subplot(2,2,3)
   plotCellData(G, sqrt(S.C'*abs(solution.cellFlux)) );
   title('Sqrt - Flux intensity')
   set_view()

subplot(2,2,4)
   h = plotCellData(G, solution.s(:,1));
   set_view()

tf = sum(porvol)/2;
dt = tf/50;
t  = 0;

while t < tf,

    solution = twophaseUpwBE(solution, dt, sources, gm, porvol, fluid, ...
                             'verbose', verbose,  'nltol', 1e-8);

    %solution = twophaseUpwFE(solution_exp, timestep, sources, gm, ...
    %                         porvol, fluid, 'verbose', verbose, 'dt', dt);     
    %                          
    %solution.s(:,1) = reorder_incompressible_twophase(solution.s(:,1), dt, sources, ...
    %                                              gm, porvol, 'transport/reorder.txt');
    %[solution.s(:,1), report] = twophaseUpwReorder(solution.s(:,1), dt, sources, gm, porvol);
    %    report
    delete(h);
    h = plotCellData(G, solution.s(:,1));
    title('Water saturation')
    drawnow

    t = t + dt;
end

function set_view()
view(3), camproj perspective, axis tight equal, camlight headlight
