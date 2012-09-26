% simpleBlackOil -- Demonstrates the Black Oil pressure solver.

clear

%--------------------------------------------------------------------------
% Reservoir description ---------------------------------------------------
%
cellDims = [ 40, 40, 10 ];
physDims = cellDims .* [100, 100, 10];

% Initial reservoir pressure
initP    =  335 * barsa();
DT       = 1000 * day()  ;

%--------------------------------------------------------------------------
% Fluid description -------------------------------------------------------
%
%pvtfile = 'pvt.txt';
%pvtfile = 'compressible.txt';
%pvtfile = 'incompressible.txt';
%pvtfile = 'idealgas.txt';
%pvtfile = 'immiscible.txt';
pvtfile = 'miscibleoil.txt';

%--------------------------------------------------------------------------
% Well controls -----------------------------------------------------------
%
useWells = true ;
useBC    = false;
useGrav  = false;
initcomp = 'Oil';
inject   = 'Water';
pHigh    = 800.0*barsa();
pLow     = 200.0*barsa();

G = cartGrid(cellDims, physDims);
G.nodes.coords(:,3) = 1950 + G.nodes.coords(:,3) +      ...
                      5*sin(.005*G.nodes.coords(:,1)) + ...
                      .02*G.nodes.coords(:,2);
G = computeGeometry(G);

diag_k    = [500, 500, 50] .* milli*darcy;
rock.perm = repmat(diag_k, [G.cells.num, 1]);
rock.poro = repmat(0.3   , [G.cells.num, 1]);

S         = computeMimeticIP(G, rock, 'Verbose', false);

%% Wells
W   = struct([]);
bc  = [];
src = [];
switch inject,
   case 'Water', inj_comp = [1, 0, 0];
   case 'Oil',   inj_comp = [0, 1, 0];
   case 'Gas',   inj_comp = [0, 0, 1];
   otherwise,
      error('Unknown injection phase: ''%s''.', inject);
end
if useWells,
   W = verticalWell(W, G, rock, 8, 8, 1:G.cartDims(3), ...
                    'Type', 'bhp', 'Val', pHigh,       ...
                    'Radius', 0.125, 'Comp_i', inj_comp);
   W = verticalWell(W, G, rock, 33, 33, ...G.cartDims(1), G.cartDims(2),      ...
                    1 : G.cartDims(3), 'Type', 'bhp', 'Val', pLow, ...
                    'Radius', 0.125, 'Comp_i', inj_comp);
   W = assembleWellSystem(G, W);
end
if useBC,
   bc = pside(bc, G, 'LEFT' , 1:cellDims(2), 1:cellDims(3), ...
              500*barsa(), 'sat', inj_comp);
   bc = pside(bc, G, 'RIGHT', 1:cellDims(2), 1:cellDims(3), ...
              400*barsa(), 'sat', inj_comp);
end

%% Define and solve pressure system
% Phase equilibrium: scale z such that phase volumes add to 1. For a given
% pressure and composition, u = inv(B(p,z))*R(p,z)'*z; and 
% alpha*u = inv(B(p,z))*R'(p,z)(alpha*z).
if true,
   fluidfile = fullfile(ROOTDIR, 'params', 'fluid', 'Data', pvtfile); 
   fluid = initBlackoilFluid(readpvt(fluidfile), 'verbose', true);
   rho_s = fluid.surfaceDensity;
   m0(3) = 0.1457;
   m0(2) = 1 - m0(3);
   z0    = bsxfun(@rdivide, m0, rho_s);
   [u, u, u, u] = fluid.pvt(initP, z0);
   alpha = sum(u, 2);
   z0    = bsxfun(@rdivide, z0, alpha);
   s0    = bsxfun(@rdivide, u,  alpha);
elseif false,
   fluid = initCompressibleFluid(500, initP, 1/initP, ...
                                 1*centi*poise());
   s0 = 1;
   z0 = 1;
else 
   fluid = initTwoCompressibleFluidsC();
   s0 = [0.5, 0.5];
   z0 = s0;
end
resSol  = initResSol (G, initP, s0, z0);
wellSol = initWellSol(W, initP);

gravity reset, gravity(useGrav);
porvol  = poreVolume(G, rock);
vd      = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);

[u, u, u, u] = fluid.pvt(initP, z0);

% Measure initial volume discrepancy term (== 0 by definition).
S.RHS.volume_discrepancy = vd(DT, u);

solver_press = @(xr, xw, p0, dt)               ...
   solveBlackOilWellSystem(xr, xw, G, rock, S, ...
                           fluid, p0, dt, 'bc', bc, 'wells', W);

% Solve pressure equation.
%
solve_p = @(xr, xw) solver_press(xr, xw, resSol.cellPressure, DT);
[resSol, wellSol] = succSubst(resSol, wellSol, solve_p,         ...
                              'Tol', 1.0e-8/day(), 'MaxIt', 40, ...
                              'Verbose', true);
[u, u, u, u] = fluid.pvt(resSol.cellPressure, resSol.z);
resSol.s = bsxfun(@rdivide, u, sum(u, 2));

%% Plot results
ax = newplot();
h  = plotCellData(G, convertTo(resSol.cellPressure, barsa()), ...
                  'EdgeColor', 'k');
view(3), camlight, grid on
title('Pressure')

[htop, htext, hs] = plotWell(G, W, 'height', 10, 'color', 'g');
colorbar
caxis(convertTo([min(resSol.cellPressure), ...
                 max(resSol.cellPressure)], barsa()));

%% Example last changed:
% $Id: simpleBlackOil.m 2343 2009-06-07 17:24:44Z bska $
