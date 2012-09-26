%% Two-Point Flux Approximation Solvers
% Use the two-point flux approximation (TPFA) method to solve the
% single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a Cartesian grid with isotropic, homogeneous permeability and two
% wells. The problem solved in this example is the same as in the tutorial
% <simpleWellExample.html "Using Peacemann Well Models">, and more a
% detailed description of the model setup is shown there.
%
% The main idea of the TPFA method is to approximate the flux v over a face
% by the difference of the cell centered pressures in the neighboring cells
% weigthed by a face transmissibility T:
%
% $$ v_{ij} = T_{ij}(p_i \textbf{--} p_j).$$
%
% The pressure in each cell is approximated by solving a linear system
% Ap = b. When ignoring wells, sources, and bc, A and b are given by
%
% $$ a_{ik} = \left\{\begin{array}{cc}
%              \sum_j t_{ij}  & \textrm{if } i=k, \\
%              -t_{ij} &  \textrm{if }\, i\neq k,
%              \end{array} \right. \quad
%              \textrm{and} \quad b_i = \int_{i} q \, dx. $$
%
% Once the pressure is known, the flux is calculated using the expression
% given above.
%
% In this example we show two different ways of building and solving a TPFA
% system. First we show the direct/classic approach using
% transmissibilities, while in the second approach we view the TPFA-method
% as a mimetic inner product and use the mimetic framework.
%
% The grid and the wells will be the same for both approaches.
%


%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 20; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G, 'Verbose', true);

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
rock.perm = repmat(100, [G.cells.num, 1]) .* 1e-3*darcy();
fluid     = initSingleFluid();
gravity off

%% Introduce wells
% We will include two wells, one rate-controlled vertical well and one
% horizontal well controlled by bottom-hole pressure. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. We need to specify ('InnerProduct', 'ip_tpf') to get the
% correct well model for TPFA.
%
% The first well is vertical well (vertical is default):
%
% * completion in cells: cellsWell1
% * controlled by production rate = 1.0  [m^3/d]
% * radius = 0.1.                        [m]
%
cellsWell1 =  1 : nx*ny : nx*ny*nz;
radius     = .1;
W = addWell(G, rock, [], cellsWell1,          ...
            'Type', 'rate', 'Val', 1.0/day(), ...
            'Radius', radius, 'InnerProduct', 'ip_tpf');

%%
%  The second well is horizontal in the 'y' direction:
%
% * completion in cells: cellsWell2
% * controlled by bottom hole pressure, bhp = 1e5 [Pa]
% * radius = 0.1                                  [m]
%
cellsWell2 = nx : ny : nx*ny;
W = addWell(G, rock,  W, cellsWell2,      ...
            'Type', 'bhp' , 'Val', 1*barsa(), ...
            'Radius', radius, 'Dir', 'y', 'InnerProduct', 'ip_tpf');


%% APPROACH 1: Direct/Classic TPFA
% Initialize solution structure with reservoir pressure equal 0. Compute
% one-sided transmissibilities for each face of the grid from input grid
% and rock properties. The harmonic averages of ones-sided
% transmissibilities are computed in the solver incompTPFA.
resSol1 = initResSol(G, 0.0);
T = computeTrans(G, rock, 'Verbose', true);

%%
% Initialize well solution structure (with correct bhp).
% No need to assemble well system (wells are added to the linear system
% inside the incompTPFA-solver).
wellSol1 = initWellSol(W, 1.0*barsa());

% Solve linear system construced from T and W to obtain solution for flow
% and pressure in the reservoir and the wells. Notice that the TPFA solver
% is different from the one used for mimetic systems.
[resSol1, wellSol1] = incompTPFA(resSol1, wellSol1,  G, T, ...
                                 fluid, 'wells', W, 'Verbose', true);


%% APPROACH 2: Mimetic with TPFA-inner product
% Initialize solution structure with reservoir pressure equal 0. Compute
% the mimetic inner product from input grid and rock properties.
resSol2 = initResSol(G, 0.0);
IP = computeMimeticIP(G, rock, 'Verbose', true, ...
                      'InnerProduct', 'ip_tpf');

%%
% Generate the components of the mimetic linear system corresponding to the
% two wells and initialize the solution structure (with correct bhp)
wellSol2 = initWellSol(W, 1.0*barsa());

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
subplot(2,2,1)
   plotCellData(G, resSol1.cellPressure ./ barsa());
   title('Pressure: direct TPFA')
   view(3), camproj perspective, axis tight off, camlight headlight
   colorbar; cax = caxis;

subplot(2,2,2)
   plotCellData(G, resSol2.cellPressure ./ barsa());
   title('Pressure: mimetic TPFA')
   view(3), camproj perspective, axis tight off, camlight headlight
   colorbar; caxis(cax);

subplot(2,2,3:4)
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
%   $Id: simpleTPFA.m 2056 2009-04-21 12:22:00Z jrn $
% </font>
% </html>
displayEndOfDemoMessage(mfilename)

