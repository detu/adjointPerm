%% Setup of iterative solution procedure
useCorrection       = true;
singleCorrection    = false;
overlap             = 5;
Nb                  = 10;
heterogen           = true;
tol                 = 1e-4; 
numSteps            = 100;
maxIt               = 50;


%% Reservoir and fluid properties
Nc       = 100;
cellDims = [ Nc, 1, 1 ];
physDims = [ 1000, 1, 1 ];  % m
G        = computeGeometry(cartGrid(cellDims, physDims));

if ~heterogen,
   diag_k = [100, 100, 100] .* milli*darcy;
   rock.perm = repmat(diag_k, [G.cells.num, 1]);
   rock.poro = repmat(0.3   , [G.cells.num, 1]);
   rock.poro(1) = 0.1;
   rock.perm(1) = 10*milli*darcy;
else
   nc = G.cells.num;
   smoother = spdiags(repmat([1/2, -1, 1/2], [nc, 1]), ...
                        -1 : 1, nc, nc);
   smoother(1,1:2) = [-1, 0];
   smoother(end,end-1:end) = [0, -1];

   poro = expm(smoother) * rand([nc, 1]);
   perm = exp(7 * poro + 1);

   rock.perm = convertFrom(perm, milli*darcy);
   rock.poro = poro;
   %load hetrock
end

%% Generate coarse grid
coarseDim = [ Nb, 1, 1 ];
p  = partitionUI(G, coarseDim);
CG = generateCoarseGrid(G, p);

%% Wells and external influence
gravity reset, gravity(false);
[W, WC, src, bc] = deal([]);
pHi     = 800.0 * barsa();
pLo     = 200.0 * barsa();
initP   = pLo;
bc      = pside(bc, G, 'LEFT', 1:1, 1:1,  pHi, 'sat', 1);
bc      = pside(bc, G, 'RIGHT', 1:1, 1:1, pLo, 'sat', 1);

%% Fluid model
% Use a three-phase model consisting of three ideal gases with identical
% properties. Only the first gas is present initially.
m0       = 1;
fluid    = Fluid(IdealGases(1*centi*poise(), 1e-5, 1*barsa));

% Phase equilibrium: scale z such that phase volumes add to 1. For a given
% pressure and composition, u = inv(B(p,z))*R(p,z)'*z; and
% alpha*u = inv(B(p,z))*R'(p,z)(alpha*z).
z0    = bsxfun(@rdivide, m0, fluid.surfaceDensity);
[u,u,u,u] = fluid.pvt(initP, z0);
alpha = sum(u, 2);
z0    = bsxfun(@rdivide, z0, alpha);
s0    = bsxfun(@rdivide, u,  alpha);
      
% define macro to compute volume discrepancy
porvol = poreVolume(G, rock);
vd     = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);


%% Finish model initialisation

% Initialize solution objects (fine grid and multiscale)
xrRef = initResSol (G, initP, s0, z0);
xrMs  = xrRef;
xwRef = initWellSol(W, initP);
xwMs  = xwRef;

% Compute time (tau) until steady state and pressure time step.
%
[mu, mu, mu, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
tau = mu(1,end) .* rock.poro(1) * physDims(1)^2 / ...
      (rock.perm(1,1) * abs(pHi - pLo));
DT  = tau / 400;

% Construct basic pressure system components.
%
s   = bsxfun(@rdivide, u, sum(u,2) );
mob = sum( bsxfun(@rdivide, fluid.relperm(s), mu), 2);
S   = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');
CS  = generateCoarseSystem(G, rock, S, CG, mob, 'bc', bc, 'src', src, 'BasisWeighting', 'poros');
W   = assembleWellSystem(G, W);
%WC  = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, 'src', src);


%% Set aside storage for all time steps
cellno        = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
itNo          = zeros(numSteps, Nb+1);
pressRef      = zeros([numSteps+1, G.cells.num]);
press         = zeros([numSteps+1, G.cells.num]);
velRef        = zeros([numSteps+1, numel(cellno)]);
vel           = zeros([numSteps+1, numel(cellno)]);
ep            = zeros(numSteps, 1);
emRef         = zeros(numSteps, 1);
emMs          = zeros(numSteps, 1);
its           = zeros(numSteps, 1);
pressRef(1,:) = xrRef.cellPressure;
press   (1,:) = xrRef.cellPressure;
velRef  (1,:) = xrRef.cellFlux;
vel     (1,:) = xrRef.cellFlux;
vdRef         = vd(DT, u);
vdMs          = vd(DT, u);


trans = computeTrans(G, rock);
%% Main loop
T = 0;
fh1 = figure(1);
fh2 = figure(2);
solver_sat   = @(xr, xw, tf, W) ...
   explicitBlackOil(G, xr, xw, tf, porvol, fluid, 'bc', bc);
for step = 1 : numSteps
   disp(['------------------------- STEP: ' num2str(step) ' ---------------------'])

   %% Compute fine-scale pressure
   S.RHS.volume_discrepancy = vdRef;
   [xrRef, xwRef] = solveMixedBO(xrRef, xwRef, G, rock, S, fluid, DT, ...
                                'bc', bc, 'tol', tol);
   xrRef = computeFacePressure(xrRef, G, trans, fluid);
   %% Compute multiscale pressure
   %S.RHS.volume_discrepancy = vdMs;
   [xrMs, xwMs, itNo(step,:)] = solveMixedBOITMS(xrMs, xwMs, G, CG, p, rock, S, CS, fluid, ...
                                   DT, 'bc', bc, 'overlap', overlap, ...
                                   'useCorrection', useCorrection, ...
                                   'tol', tol, 'singleCorrection', singleCorrection);
   xrMs = computeFacePressure(xrMs, G, trans, fluid);
 
   %% Plot pressures
   FontSz = 8;
   set(0,'CurrentFigure',fh1);
   subplot(3,2,[1 3]),
      X = G.cells.centroids(:,1);
      plot(X, convertTo([xrRef.cellPressure, xrMs.cellPressure], barsa), ...
          'LineWidth', 2)
      v = axis; axis([v(1:2), convertTo(min(pLo,pHi),barsa)-1, convertTo(max(pHi,pLo), barsa)])
      title('Cell Pressure', 'FontSize', FontSz)
      legend('Reference', 'MsMFEM'), set(gca, 'FontSize', FontSz)
  
 
   %% Compute fluid transport
   xrRef = solver_sat(xrRef, xwRef, DT, W );
   xrMs  = solver_sat(xrMs , xwMs , DT, WC);
 
   % Compute volume discrepancies
   [u, u, u, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
   vdRef        = vd(DT, u);
   [u, u, u, u] = fluid.pvt(xrMs.cellPressure, xrMs.z);
   vdMs         = vd(DT, u);

   massRef = bsxfun(@times, xrRef.z, fluid.surfaceDensity);
   massMs  = bsxfun(@times, xrMs .z, fluid.surfaceDensity);
  
   relPDiff = (xrMs.cellPressure - xrRef.cellPressure)/(max(xrRef.cellPressure) - min(xrRef.cellPressure));
   ep(step) = norm(relPDiff, inf);

   %% Various plots
   FontSz = 8;
   subplot(3,2,5),
      %fp = convertTo(relPDiff, barsa);
      plot(X,relPDiff, 'LineWidth', 2)
      title('Relative difference in cell pressure')
      set(gca, 'FontSize', FontSz);
  
   subplot(3,2,6),
      relMassErrorRef = DT * (vdRef./porvol);
      relMassErrorMs  = DT * (vdMs./porvol);
   
      emRef(step)  = sum(abs(relMassErrorRef))/Nc;
      emMs(step)   = sum(abs(relMassErrorMs))/Nc;
   
   cla, hold on
      plot(X, (massRef-massMs)./massRef, 'LineWidth', 2);
      title('Relative discrepancy in mass');
      %set(gca, 'XTick', []),
      set(gca, 'FontSize', FontSz);

   subplot(3,2,2),
      ff = convertTo([xrRef.faceFlux, xrMs.faceFlux], meter^3/day);
      %plot(X,ff(1 : Nc, :)*[1 -1]'/norm(ff(1 : Nc, 1)), 'LineWidth', 2)
      plot(X,ff(1 : Nc, :), 'LineWidth', 2)
      title('Total flux'), 
      %set(gca, 'XTick', []),
      set(gca, 'FontSize', FontSz);


   subplot(3,2,4),
   cla, hold on
      plot(X, massRef, 'LineWidth', 2);
      plot(X, massMs, '--','LineWidth', 2);
      title('Total mass');
      set(gca, 'FontSize', FontSz);

   set(0,'CurrentFigure',fh2);
      bar(itNo(1:step,:),'stacked');
      set(gca,'XLim',[1 max(step,2)]);
      
   T = T + DT;

   %% Store this time step
   drawnow
   pressRef(step+1,:) = xrRef.cellPressure;
   velRef  (step+1,:) = xrRef.cellFlux;
   press   (step+1,:) = xrMs .cellPressure;
   vel     (step+1,:) = xrMs .cellFlux;
end
