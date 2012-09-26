function [resSol, wellSol] = blackOilQuarterFiveSpot
% blackOilQuarterFiveSpot -- Solve synthetic, quarter-five-spot 
%                            reservoir using black oil model
%
% SYNOPSIS:
%   [resSol, wellSol] = blackOilQuarterFiveSpot
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   resSol  - Reservoir solution structure at steady state.
%
%   wellSol - Well solution structure at steady state.
%
% SEE ALSO:
%   computeMimeticIP, assembleWellSystem, solveBlackOilWellSystem.

% $Id:$

c = onCleanup(@()exportWorkspace); %#ok

% Whether or not to emit convergence history during pressure solves
verbose  = true;  %#ok

%--------------------------------------------------------------------------
% Reservoir description ---------------------------------------------------
%
cellDims = [ 40, 40, 1 ];
physDims = cellDims      .* [10, 10, 10];

useGrav = false;

% Initial reservoir pressure
initP    = 420e5;

%--------------------------------------------------------------------------
% Fluid description -------------------------------------------------------
%
%pvtfile = 'pvt.txt';
pvtfile = 'compressible.txt';
%pvtfile = 'incompressible.txt';
%pvtfile = 'idealgas.txt';

%--------------------------------------------------------------------------
% Well controls -----------------------------------------------------------
%
useWells   = true;
initcomp   = 'Oil';
inject     = 'Water';
injRate    = 1000.0 / 86400;   % [m^3/s]
prodPress  =  420.0e5;         % Pa

%--------------------------------------------------------------------------
% Controls on numerical methods -------------------------------------------
%
DT   = 1.0 * 86400;          % Pressure time step restriction [s]

%--------------------------------------------------------------------------
%% Reservoir and fluid properties -----------------------------------------
%
PVTTAB     = readpvt(fullfile(ROOTDIR, 'params', 'fluid', 'Data', pvtfile));

G          = cartGrid(cellDims, physDims);
G          = computeGeometry(G);

rock.perm  = repmat([500.0, 500.0, 50.0], [G.cells.num, 1]);   % mD
rock.poro = repmat(                0.3 , [G.cells.num, 1]);

rock.perm  = rock.perm .* darcy() / 1000;  % mD -> m^2

S          = computeMimeticIP(G, rock, 'Verbose', false);

%--------------------------------------------------------------------------
%% Wells and external influence -------------------------------------------
%
W  = struct([]);
bc = [];
switch inject,
   case 'Water', inj_comp = [1, 0, 0];
   case 'Oil',   inj_comp = [0, 1, 0];
   case 'Gas',   inj_comp = [0, 0, 1];
   otherwise,
      error('Unknown injection phase: ''%s''.', inject);
end
if useWells,
   W = addWell(G, rock, W, 1,      ...
               'Type',   'rate',   ...
               'Val',    injRate,  ...
               'Radius', 0.10,     ...
               'Comp_i', inj_comp);
   W = addWell(G, rock, W, G.cells.num, ...
               'Type',   'bhp',         ...
               'Val',    prodPress,     ...
               'Radius', 0.10);
   W = assembleWellSystem(G, W);
else
   bc = fluxside(bc, G, 'LEFT',  1:cellDims(2), 1:cellDims(3), ...
                 injRate, 'sat', inj_comp);
   bc = pside   (bc, G, 'RIGHT', 1:cellDims(2), 1:cellDims(3), ...
                 prodPress, 'sat', inj_comp);
end

%--------------------------------------------------------------------------
%% Finish model initialisation --------------------------------------------
%
resSol  = initResSol (G, initP, [0, 0, 0]);
wellSol = initWellSol(W, initP);

switch initcomp,
   case 'Water', resSol.z(:,1) = 1;
   case 'Oil',   resSol.z(:,2) = 1;
   case 'Gas',   resSol.z(:,3) = 1;
   otherwise,
      error('Unknown initial reservoir phase: ''%s''.', initcomp);
end

gravity(useGrav);

[fluid, resSol] = pvt(PVTTAB, resSol, gravity());

% Set correct mass values following the new saturations in 'resSol'.
%
resSol.z = resSol.s ./ fluid.B;
porvol   = poreVolume(G, rock);

vd = @(x, dt, fluid) ((sum(x.z .* fluid.B, 2) - 1) .* porvol ./ dt);

% Measure initial volume discrepancy term (== 0 by definition).
%
solver_press = @(xr, xw, p0, dt, S)                  ...
   solveBlackOilWellSystem(xr, xw, G, rock, S, W, ...
                           PVTTAB, p0, dt, 'bc', bc);
solver_sat   = @(xr, xw, tf, fluid) ...
   implicitBlackOil(G, xr, xw, tf, porvol,     ...
                    fluid.B, fluid.mu, PVTTAB, ...
                    'bc', bc, 'wells', W);

plot_pd = @(x)   plot(1 : G.cartDims(1), ...
                      x.cellPressure(1 : G.cartDims(1) + 1 : end) ./ barsa());
plot_sd = @(x,c) plotyy(1 : G.cartDims(1),                      ...
                        x.s(1 : G.cartDims(1) + 1 : end, c(1)), ...
                        1 : G.cartDims(1),                      ...
                        x.s(1 : G.cartDims(1) + 1 : end, c(2)));

%--------------------------------------------------------------------------
%% Run model until convergence --------------------------------------------
%
T = 0;
for k = 1 : 400,
   % Volume discrepancy source attempt to correct the pressue such that the
   % the sum of phase volume fractions will be one at the next time step.
   %
   S.RHS.volume_discrepancy = vd(resSol, DT, fluid);

   % Solve pressure equation
   %
   solve_p = @(xr, xw) solver_press(xr, xw, resSol.cellPressure, DT, S);
   [resSol, wellSol] = succSubst(resSol, wellSol, solve_p,   ...
                                 'Tol', 1.0e-8/86400, 'MaxIt', 40, ...
                                 'Verbose', true);
   subplot(3,1,1),
      plot_pd(resSol);
      title('Pressure -- diagonal [bar]');

   % Solve saturation equation
   %
   [fluid, resSol] = pvt(PVTTAB, resSol, gravity());
   resSol          = solver_sat(resSol, wellSol, DT, fluid);
   resSol.s(:)     = resSol.z .* fluid.B;

   subplot(3,1,2),
      ax = plot_sd(resSol, [1, 2]);
      title('Densities -- diagonal');
      axis(ax(1), [1, G.cartDims(1), 0, 1]);
      axis(ax(2), [1, G.cartDims(1), 0, 1]);

   subplot(3,1,3),
      surf(reshape(resSol.s(:,1), G.cartDims(1:2))); view(45, 55);
      title('Saturation -- reservoir');

   drawnow

   T = T + DT;
end
