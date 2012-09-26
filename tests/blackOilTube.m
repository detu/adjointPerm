function [resSol, wellSol] = blackOilTube
% blackOilTube -- Solve synthetic, tubular reservoir using black oil model
%
% SYNOPSIS:
%   [resSol, wellSol] = blackOilTube
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

% $Id: blackOilTube.m 2338 2009-06-05 17:19:30Z bska $

cleanup = onCleanup(@()exportWorkspace); %#ok

% Whether or not to emit convergence history during pressure solves
verbose  = true;  %#ok

%--------------------------------------------------------------------------
% Reservoir description ---------------------------------------------------
%
cellDims = [ 100, 1, 1 ];
physDims = cellDims        .* [1, 1, 1];

% Initial reservoir pressure [Pa]
initP    = 100*barsa();

%--------------------------------------------------------------------------
% Fluid description -------------------------------------------------------
%
%pvtfile = 'pvt.txt';
%pvtfile = 'compressible.txt';
pvtfile = 'incompressible.txt';
%pvtfile = 'idealgas.txt';
%pvtfile = 'immiscible.txt';

%--------------------------------------------------------------------------
% Well controls -----------------------------------------------------------
%
useWells   = true;
%initcomp   = 'Oil';
inject     = 'Water';
injPress   = 200.0*barsa();  % [Pa]
prodPress  = 100.0*barsa();  % [Pa]

%--------------------------------------------------------------------------
% Controls on numerical methods -------------------------------------------
%
TF   = 10*day();
DT   = day() / 10000;     % Pressure time step restriction [s]

%--------------------------------------------------------------------------
%% Reservoir and fluid properties -----------------------------------------
%

G          = cartGrid(cellDims, physDims);
G          = computeGeometry(G);


diag_k     = [500.0, 500.0, 50.0] .* (darcy() / 1000);  % mD -> m^2
rock.perm  = repmat(diag_k, [G.cells.num, 1]);
rock.poro = repmat(0.3 ,   [G.cells.num, 1]);

S          = computeMimeticIP(G, rock, 'Verbose', false);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Finish model initialisation --------------------------------------------
%
model = 'eclipse';
%model = 'JuanesPatzek';

switch model
   case 'eclipse'
      fluidfile = fullfile(ROOTDIR, 'params', 'fluid', 'Data', pvtfile);
      fluid = initBlackoilFluid(readpvt(fluidfile), 'verbose', true);
      rho_s = fluid.surfaceDensity;
      m0(3) = 0.5;
      m0(1) = 0.5;
      m0(2) = 1-m0(3)-m0(1);
      z0    = bsxfun(@rdivide, m0, rho_s);
      [u, u, u, u] = fluid.pvt(initP, z0);
      alpha = sum(u, 2);
      z0    = bsxfun(@rdivide, z0, alpha);
      s0    = bsxfun(@rdivide, u,  alpha);
      
      inj_comp    = z0;
      
   case 'singlefluid'
      fluid = initCompressibleFluid(500, initP, 1/initP, 1*centi*poise());
      s0 = 1;
      z0 = 1;
   case 'twofluids'
      fluid = initTwoCompressibleFluidsC();
      s0 = [0.5, 0.5];
      z0 = s0;
   case 'JuanesPatzek'
      
      fluid = initSimpleThreephaseCompressibleFluid([1000, 700, 70], ...
                                                    100*barsa(), ...
                                                    [0,0,2/(100*barsa())], ...
                                                    [0.875, 2.0, 0.03]*centi*poise());
      rho_s = fluid.surfaceDensity;
      m0(3) = 0.5;
      m0(1) = 0.5;
      m0(2) = 1-m0(3)-m0(1);
      z0    = bsxfun(@rdivide, m0, rho_s);
      [u, u, u, u] = fluid.pvt(initP, z0);
      alpha = sum(u, 2);
      z0    = bsxfun(@rdivide, z0, alpha);
      s0    = bsxfun(@rdivide, u,  alpha);
      
      inj_comp    = z0;
      
   otherwise
      error('Huh!? %s --- What fluid model is that?', model);
end

%% Wells and external influence -------------------------------------------
%
W = struct([]);
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
               'Type',   'bhp',   ...
               'Val',    injPress,  ...
               'Radius', 0.10,     ...
               'Comp_i', inj_comp);
   W = addWell(G, rock, W, G.cells.num, ...
               'Type',   'bhp',         ...
               'Val',    prodPress,     ...
               'Radius', 0.10);
   W = assembleWellSystem(G, W);
else
   bc = fluxside(bc, G, 'LEFT', 1:cellDims(2), 1:cellDims(3),  ...
                 injRate, 'sat', inj_comp);
   bc = pside   (bc, G, 'RIGHT', 1:cellDims(2), 1:cellDims(3), ...
                 prodPress, 'sat', inj_comp);
end






resSol  = initResSol (G, initP, s0, z0);
wellSol = initWellSol(W, initP);

gravity(false);
porvol   = poreVolume(G, rock);
vd       = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);


[u, u, u, u] = fluid.pvt(initP, z0);
% resSol.s = bsxfun(@rdivide, u, sum(u, 2));

% Measure initial volume discrepancy term (== 0 by definition).
S.RHS.volume_discrepancy = vd(DT, u);


% Measure initial volume discrepancy term (== 0 by definition).
%
solver_press = @(xr, xw, p0, dt, S)                  ...
   solveBlackOilWellSystem(xr, xw, G, rock, S, ...
                           fluid, p0, dt, 'bc', bc, 'wells', W);
solver_sat   = @(xr, xw, tf) ...
   explicitBlackOil(G, xr, xw, tf, porvol,     ...
                    fluid, ...
                    'bc', bc, 'wells', W);

%--------------------------------------------------------------------------
%% Run model until convergence --------------------------------------------
%
T = 0;
for k = 1 : TF/DT,

   % Solve pressure equation
   %
   solve_p = @(xr, xw) solver_press(xr, xw, resSol.cellPressure, DT, S);
   [resSol, wellSol] = succSubst(resSol, wellSol, solve_p,   ...
                                 'Tol', 1.0e-8, 'MaxIt', 40, ...
                                 'Verbose', false);

   subplot(3,2,1), plot(resSol.cellPressure./barsa()); title('Cell pressure [bar]');

   % Solve saturation equation
   %
   %[fluid, resSol] = pvt(PVTTAB, resSol, gravity());
   resSol          = solver_sat(resSol, wellSol, DT);
   %resSol.s(:)     = resSol.z .* fluid.B;

   %resSol.s(:) = bsxfun(@rdivide, u, sum(u,2));
   
   [c, rho, mu, u, R] = fluid.pvt(resSol.cellPressure, resSol.z);%#ok
   rhoS = fluid.surfaceDensity;
   mass = bsxfun(@times, resSol.z, rhoS);
   subplot(3,2,2), plot(mass(:,1)); title('Water [kg]');
   subplot(3,2,4), plot(mass(:,2)); title('Oil [kg]');
   subplot(3,2,6), plot(mass(:,3)); title('Gas [kg]');

   % Volume discrepancy source attempt to correct the pressue such that the
   % the sum of phase volume fractions will be one at the next time step.
   %
   S.RHS.volume_discrepancy = vd(DT, u);
   subplot(3,2,3), plot(DT .* S.RHS.volume_discrepancy); title('Volume error');
   
   drawnow

   T = T + DT;
end

