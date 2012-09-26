function [resSol, wellSol] = blackOilTubeTPF
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

% $Id: blackOilTube.m 668 2008-09-05 12:27:32Z bska $

c = onCleanup(@()exportWorkspace); %#ok

% Whether or not to emit convergence history during pressure solves
verbose  = true;  %#ok

%--------------------------------------------------------------------------
% Reservoir description ---------------------------------------------------
%
cellDims = [ 1000, 1, 1 ];
physDims = cellDims        .* [1, 1, 1];

% Initial reservoir pressure
initP    = 100;

%--------------------------------------------------------------------------
% Fluid description -------------------------------------------------------
%
%pvtfile = 'pvt.txt';
%pvtfile = 'compressible.txt';
%pvtfile = 'incompressible.txt';
pvtfile = 'idealgas.txt';

%--------------------------------------------------------------------------
% Well controls -----------------------------------------------------------
%
useWells   = true;
initcomp   = 'Gas';
inject     = 'Gas';
CFL        = 0.01;
injRate    = 300.0;
prodPress  = 100.0;

%--------------------------------------------------------------------------
% Controls on numerical methods -------------------------------------------
%
pTol = 5.0e-7;      % Absolute convergence tolerance, pressure
vTol = 5.0e-6;      % Absolute convergence tolerance, flux
zTol = 5.0e-9;      % Absolute convergence tolerance, mass distribution

DT   = 1 / 100;     % Pressure time step restriction

%--------------------------------------------------------------------------
%% Reservoir and fluid properties -----------------------------------------
%
PVTTAB     = readpvt(fullfile(SAMSIMROOTDIR, 'fluid', 'Data', pvtfile));

G          = cartGrid(cellDims, physDims);
G          = computeGeometry(G);


rock.perm  = repmat([500.0, 500.0, 50.0], [G.cells.num, 1]);   % mD
rock.poro = repmat(                0.3 , [G.cells.num, 1]);

S          = computeMimeticIP(G, rock,          ...
                                   'InnerProduct', 'ip_tpf', ...
                                   'Verbose', true, ...
                                   'Gravity', false);

%--------------------------------------------------------------------------
%% Wells and external influence -------------------------------------------
%
W = struct([]);
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
   S = fluxside(S, G, 'LEFT',  1:1, 1:1, injRate  );
   S = pside   (S, G, 'RIGHT', 1:1, 1:1, prodPress);
end

%--------------------------------------------------------------------------
%% Finish model initialisation --------------------------------------------
%
resSol  = initResSol (G, initP);
wellSol = initWellSol(W, initP);

switch initcomp,
   case 'Water', resSol.z(:,1) = 1;
   case 'Oil',   resSol.z(:,2) = 1;
   case 'Gas',   resSol.z(:,3) = 1;
   otherwise,
      error('Unknown initial reservoir phase: ''%s''.', initcomp);
end

[fluid, resSol]  = pvt(PVTTAB, resSol, S.constants.gravity);

% Set correct mass values following the new saturations in 'resSol'.
%
resSol.z         = resSol.s ./ fluid.B;
porvol           = poreVolume(G, rock);

% Measure initial volume discrepancy term (== 0 by definition).
%
S.RHS.volume_discrepancy = (sum(resSol.z .* fluid.B, 2) - 1) .* porvol ./ DT;

%--------------------------------------------------------------------------
%% Run model until convergence --------------------------------------------
%
T  = 0;
converged = false;
while ~converged,

   % Solve pressure equation
   %
   resV = 1;

   p0 = resSol.cellPressure;
   v0 = resSol.cellFlux;
   z0 = resSol.z;

   while resV > 1.0e-8,
      cf0 = resSol.cellFlux;
      %********************************************************************
      [resSol, wellSol, S] = solveTPFBlackOilWellSystem      ...
         (resSol, wellSol, G, rock, S, W, PVTTAB, p0, DT, ...
          'Verbose', false);
      %********************************************************************
      resV = norm(resSol.cellFlux - cf0, inf);
      subplot(3,1,1)
      plot(resSol.cellPressure) %, axis([1, G.cells.num, 50, 1200]);
      drawnow

      dispif(verbose, 'Residual in velocity: %5.5e\n', resV);
   end
   [fluid, resSol] = pvt(PVTTAB, resSol, S.constants.gravity);

   % Solve saturation equation
   %
   [gm, q]  = inflow_bo(G, S, W, resSol, wellSol, PVTTAB);
   resSol.z = blackoilUpwFE(resSol.z, DT, q, gm, porvol, ...
                            fluid, PVTTAB, 'CFL', CFL);

   subplot(3,1,2), plotyy(1:G.cells.num, resSol.z(:,1), ...
                          1:G.cells.num, resSol.z(:,3))

   % Volume discrepancy source attempt to correct the pressue such that the
   % the sum of phase volume fractions will be one at the next time step.
   %
   S.RHS.volume_discrepancy = (sum(resSol.z .* fluid.B, 2) - 1) ./ DT .* porvol;
   subplot(3,1,3), plot(DT .* S.RHS.volume_discrepancy)

   % Check for convergence (i.e., steady state)
   %
   dp = norm(p0 - resSol.cellPressure)               / sqrt(G.cells.num);
   dv = sqrt(sum(S.C.'*((v0 - resSol.cellFlux).^2))) / sqrt(G.cells.num);
   dz = norm(z0 - resSol.z           )               / sqrt(G.cells.num);
   converged = (dp < pTol) & (dv < vTol) & (dz < zTol);

   T = T + DT;
end
