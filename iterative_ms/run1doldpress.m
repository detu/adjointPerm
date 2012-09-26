%% Setup of iterative solution procedure
heterogen           = false;
tol                 = 1e-12; 
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
   rock.poro(10) = 0.1;
   rock.perm(10) = 10*milli*darcy;
else
   nc = G.cells.num;
   smoother = spdiags(repmat([1/2, -1, 1/2], [nc, 1]), ...
                        -1 : 1, nc, nc);
   smoother(1,1:2) = [-1, 0];
   smoother(end,end-1:end) = [0, -1];

   poro = expm(smoother) * rand([nc, 1]);
   perm = exp(5 * poro + 1);

   rock.perm = convertFrom(perm, milli*darcy);
   rock.poro = poro;
   %load hetrock
end

%% Wells and external influence
gravity reset, gravity(false);
pHi     = 10.0 * barsa();
pLo     =  1.0 * barsa();
initP   = pLo;[W, WC, src, bc] = deal([]);
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
xrNew = xrRef;
xwRef = initWellSol(W, initP);
xwNew = xwRef;

% Compute time (tau) until steady state and pressure time step.
%
[mu, mu, mu, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
tau = mu(1,end) .* rock.poro(1) * physDims(1)^2 / ...
      (rock.perm(1,1) * (pHi - pLo));
DT  = tau / 400;

% Construct basic pressure system components.
%
s   = bsxfun(@rdivide, u, sum(u,2) );
mob = sum( bsxfun(@rdivide, fluid.relperm(s), mu), 2);
S   = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');
W   = assembleWellSystem(G, W);

%% Set aside storage for all time steps
cellno        = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
pressRef      = zeros([numSteps+1, G.cells.num]);
press         = zeros([numSteps+1, G.cells.num]);
velRef        = zeros([numSteps+1, numel(cellno)]);
vel           = zeros([numSteps+1, numel(cellno)]);
ep            = zeros(numSteps, 1);
emRef         = zeros(numSteps, 1);
emNew         = zeros(numSteps, 1);
its           = zeros(numSteps, 1);
pressRef(1,:) = xrRef.cellPressure;
press   (1,:) = xrRef.cellPressure;
velRef  (1,:) = xrRef.cellFlux;
vel     (1,:) = xrRef.cellFlux;
vdRef         = vd(DT, u);
vdNew         = vd(DT, u);

%% Main loop
T = 0;
trans = computeTrans(G, rock);
solver_sat   = @(xr, xw, tf, W) ...
   explicitBlackOil(G, xr, xw, tf, porvol, fluid, 'bc', bc);

solver_press = @(xr, xw, p0, dt, S)                  ...
   solveBlackOilWellSystem(xr, xw, G, rock, S, ...
                           fluid, p0, dt, 'bc', bc, 'wells', W);

for step = 1 : numSteps
   disp(['------------------------- STEP: ' num2str(step) ' ---------------------'])

   %% Compute fine-scale pressure, old formulation
   S.RHS.volume_discrepancy = vdRef;
   solve_p = @(xr, xw) solver_press(xr, xw, xrRef.cellPressure, DT, S);
   [xrRef, xwRef] = succSubst(xrRef, xwRef, solve_p,   ...
                                 'Tol', tol*1e-8, 'MaxIt', 40, ...
                                 'Verbose', false);

   %% Compute fine-scale pressure, residual formulation
   S.RHS.volume_discrepancy = vdNew;
   [xrNew, xwNew] = solveMixedBO(xrNew, xwNew, G, rock, S, fluid, DT, ...
                                'bc', bc, 'tol', tol);
   xrNew = computeFacePressure(xrNew, G, trans, fluid);

   %% Plot discrepancies in pressure solutions
   FontSz = 8;
   figure(1);
   subplot(3,2,1),
      Xc = G.cells.centroids(:,1);
      relcPdiff = (xrNew.cellPressure - xrRef.cellPressure)/ ...
         (max(xrRef.cellPressure) - min(xrRef.cellPressure));
      ep(step) = norm(relcPdiff, inf);
      plot(Xc,relcPdiff, 'LineWidth', 2)
      title('Discrepancy, cell pressure')
      set(gca, 'FontSize', FontSz);

   subplot(3,2,2),
      relfPdiff = (xrNew.facePressure - xrRef.facePressure)/ ...
         (max(xrRef.facePressure) - min(xrRef.facePressure));
      plot(relfPdiff, 'LineWidth', 2)
      title('Discrepancey, face pressure')
      set(gca, 'FontSize', FontSz);

   subplot(3,2,3),
      plot( (xrNew.cellFlux-xrRef.cellFlux)/max(xrRef.cellFlux) );
      title('Discrepancy, cell flux');
      set(gca, 'FontSize', FontSz);
   
   subplot(3,2,4),
      plot( (xrNew.faceFlux-xrRef.faceFlux)/max(xrRef.faceFlux) );
      title('Discrepancy, cell flux');
      set(gca, 'FontSize', FontSz);
  
   %% Compute fluid transport
   xrRef = solver_sat(xrRef, xwRef, DT, W );
   xrNew = solver_sat(xrNew, xwNew, DT, W );
 
   % Compute volume discrepancies
   [u, u, u, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
   vdRef        = vd(DT, u);
   [u, u, u, u] = fluid.pvt(xrNew.cellPressure, xrNew.z);
   vdNew        = vd(DT, u);

   massRef = bsxfun(@times, xrRef.z, fluid.surfaceDensity);
   massNew = bsxfun(@times, xrNew.z, fluid.surfaceDensity);
  

   %% Various plots
   FontSz = 8;
  
   subplot(3,2,5),
   cla, hold on
      plot(Xc, massRef, 'LineWidth', 2);
      plot(Xc, massNew, '--','LineWidth', 2);
      set(gca,'YLim',[0 1.1*max(massRef(:))]);
      title('Total mass');
      set(gca, 'FontSize', FontSz);

   subplot(3,2,6),
   cla, hold on
      relMassErrorRef = DT * (vdRef./porvol);
      relMassErrorNew = DT * (vdNew./porvol);
      emRef(step)  = sum(abs(relMassErrorRef))/Nc;
      emNew(step)   = sum(abs(relMassErrorNew))/Nc;
      plot(Xc, (massRef-massNew)./massRef, 'LineWidth', 2);
      title('Relative discrepancy in mass');
      set(gca, 'FontSize', FontSz);

   T = T + DT;

   %% Store this time step
   drawnow
   pressRef(step+1,:) = xrRef.cellPressure;
   velRef  (step+1,:) = xrRef.cellFlux;
   press   (step+1,:) = xrNew.cellPressure;
   vel     (step+1,:) = xrNew.cellFlux;
end
