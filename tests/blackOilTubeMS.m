function resSol = blackOilTubeMS
% blackOilQuarterFiveSpot
%
% implisitt løser?
% blandbar strøm
% rapportér massebalanse.
% OIP,WIP,GIP,etc
%

c = onCleanup(@exportWorkspace);  %#ok

%pvtfile    = 'pvt.txt';
%pvtfile    = 'compressible.txt';
%pvtfile    = 'incompressible.txt';
pvtfile    = 'idealgas.txt';

verbose    = true;

PVTTAB     = readpvt(fullfile(SAMSIMROOTDIR, 'fluid', 'Data', pvtfile));


cellDims   = [ 100, 1, 1];
physDims   = [ 100, 1, 1];
G          = cartGrid(cellDims, physDims);
G          = computeGeometry(G);


rock.perm  = repmat([500.0, 500.0, 50.0], [G.cells.num, 1]);   % mD
rock.poro = repmat(                0.3,  [G.cells.num, 1]);

S          = computeMimeticIP(G, rock,                  ...
                                   'Type',    'comp_hybrid', ...
                                   'Verbose', false,         ...
                                   'Gravity', false);

p  = partitionUI(G, [10, 1, 1]);
p  = processPartition  (G, p, 'Verbose', verbose);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);
CS = generateCoarseSystem(G, rock, S, CG, ...
                          'BasisWeighting', 'poros', ...
                          'Verbose', verbose);

%% Wells
W = struct([]);
W = addWell(G, rock, W, 1,    ...
            'Type',   'rate', ...
            'Val',    100.0,  ...
            'Radius', 0.10);
W = addWell(G, rock, W, G.cells.num, ...
            'Type',   'bhp',         ...
            'Val',    100,           ...
            'Radius', 0.10);
W = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, rock, W);

% Splittesteg
DT      = 1/100;

initP   = 100; %bar

resSol  = initResSol (G, initP);
wellSol = initWellSol(W, initP);

resSol.z(:,3) = 1;

pressure_prev = resSol.cellPressure;

[fluid, resSol]  = pvt(PVTTAB, resSol, S.constants.gravity);

% Set correct mass values following the new saturations in 'fluid'.
resSol.z         = resSol.s ./ fluid.B;

% Measure initial volume discrepancy term (== 0 by definition).
source = (sum(resSol.z .* fluid.B, 2) - 1) .* poreVolume(G, rock) ./ DT;
S .RHS.volume_discrepancy = source;
CS.RHS.volume_discrepancy = accumarray(p, source);

cc = accumarray(p,G.cells.centroids(:,1)) ./ accumarray(p,1);

% Splittesteg

T  = 0;
while true,
   %
   % Solve pressure
   %
   resV = 1;
   while resV > 0.5e-5,
      cf0 = resSol.cellFlux;

      [resSol, wellSol] = solveBlackOilWellSystemMS ...
         (resSol, wellSol, G, CG, rock,             ...
         S, CS, W, PVTTAB, pressure_prev, DT);

      resV = norm(resSol.cellFlux - cf0, inf);
      subplot(3,1,1)
      plot(resSol.cellPressure);
      drawnow

      dispif(verbose, 'Residual in velocity: %5.5e\n', resV);
   end
   resSol.cellPressure = interp1(cc, resSol.blockPressure, ...
      G.cells.centroids(:,1),   ...
      'spline', 'extrap');
   [fluid, resSol] = pvt(PVTTAB, resSol, S.constants.gravity);

   % Solve saturation equation
   [gm, sources, porvol] = inflow(G, rock, S, W, resSol, wellSol);
   resSol.z = blackoilUpwFE(resSol.z, DT, sources, gm, porvol, ...
                            fluid, PVTTAB);

   subplot(3,1,2),
   ax = plotyy(1:G.cells.num, resSol.z(:,1), 1:G.cells.num, resSol.z(:,3));
   axis(ax(2), [0, 100, 0, 1]);

   T = T + DT;

   pressure_prev = resSol.cellPressure;

   % Volume discrepancy source attempt to correct the pressue
   % such that the the sum of phase volume fractions will be one
   % at the next time step.
   %
   source = (sum(resSol.z .* fluid.B, 2) - 1) ./ DT .* porvol;
   S .RHS.volume_discrepancy = source;
   CS.RHS.volume_discrepancy = accumarray(p, source);
   subplot(3,1,3), plot(DT .* CS.RHS.volume_discrepancy)
end

