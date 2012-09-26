% Set variables:
useCorrection       = true;
% correctionOverlap   = 110;
% %useMSCIter          = false;
% mscIterNum          = 1;
% mscTol              = 1e-8;
numSteps            = 40;


% Whether or not to emit convergence history during pressure solves
verbose  = true;

%--------------------------------------------------------------------------
% Reservoir description ---------------------------------------------------
%
physDims = [ 1000, 1, 1 ];  % m

Nc = 100;
Nb = 5  ;
cellDims  = [ Nc, 1, 1 ];
coarseDim = [ Nb, 1, 1 ];

X = linspace(0, physDims(1), Nc);

%--------------------------------------------------------------------------
% Well controls -----------------------------------------------------------
% 
%
useWells = false;
useBC    = true ;
useSRC   = false;
useGrav  = false;
CFL      =   0.01;
pHi      =  400.00 * barsa();
pLo      =  350.00 * barsa();


initP = pLo;

%--------------------------------------------------------------------------
% Controls on numerical methods -------------------------------------------
%
%pTol = 5.0e-7;
vTol = 5.0e-8;      % Absolute convergence tolerance, flux
zTol = 5.0e-8;      % Absolute convergence tolerance, mass distribution

pTol  = 5.0e-8;     % Absolute convergence tolerance, pressure
maxIt = 20;

%--------------------------------------------------------------------------
%% Reservoir and fluid properties -----------------------------------------
%
G      = computeGeometry(cartGrid(cellDims, physDims));
cellno = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';

homogeneous = true;


if homogeneous,
   diag_k = [100, 100, 100] .* milli*darcy;
   rock.perm = repmat(diag_k, [G.cells.num, 1]);
   rock.poro = repmat(0.3   , [G.cells.num, 1]);
else
   %load hetrock
   nc = G.cells.num;
   smoother = spdiags(repmat([1/2, -1, 1/2], [nc, 1]), ...
                        -1 : 1, nc, nc);
   smoother(1,1:2) = [-1, 0];
   smoother(end,end-1:end) = [0, -1];

   poro = expm(smoother) * rand([nc, 1]);
   perm = exp(5 * poro + 1);

   rock.perm = convertFrom(perm, milli*darcy);
   rock.poro = poro;
end

p  = partitionUI(G, coarseDim);
CG = generateCoarseGrid(G, p);

%--------------------------------------------------------------------------
%% Wells and external influence -------------------------------------------
%
[W, WC, src, bc] = deal([]);
if useWells,
   W = addWell(G, rock, W, 1, ...:100:400,      ...
               'Type',   'bhp',    ...
               'Val',    pHi,      ...
               'Radius', 0.10,     ...
               'Comp_i', 1);
   W = addWell(G, rock, W, G.cells.num, ...:100:400, ...
               'Type',   'bhp',         ...
               'Val',    pLo,           ...
               'Radius', 0.10, 'Comp_i', 0);
end
if useBC,
   bc = pside(bc, G, 'LEFT', 1:1, 1:1,  pHi, 'sat', [1 0 0]);
   bc = pside(bc, G, 'RIGHT', 1:1, 1:1, pLo, 'sat', [1 0 0]);
end
if useSRC,
   src = addSource(src, 1          ,  1.0, 'sat', 1);
   src = addSource(src, G.cells.num, -1.0, 'sat', 0);
end

%% Define and solve pressure system
% Phase equilibrium: scale z such that phase volumes add to 1. For a given
% pressure and composition, u = inv(B(p,z))*R(p,z)'*z; and
% alpha*u = inv(B(p,z))*R'(p,z)(alpha*z).
inj_comp = [1, 0, 0];
m0    = [0, .05, .95];

%fluid = initSimpleThreephaseCompressibleFluid([1000, 700, 70], initP, ...
%                                                    [0, 0, 2/(initP)], ...
%                                                    [0.875, 2.0, 0.03]*centi*poise());
                                                
%fluid = ConstantCompressibilityFluids([1e-10 1e-7 0], [1000 700 1], 1*barsa);                                                
f1 = IdealGases(1*centi*poise(), 1e-5);
f2 = ConstantCompressibilityFluids(1e-8, 700 , 1*barsa);
fluid = Fluid(f1,f1,f2);

z0    = bsxfun(@rdivide, m0, fluid.surfaceDensity);
      [u, u, u, u] = fluid.pvt(initP, z0);
      alpha = sum(u, 2);
      z0    = bsxfun(@rdivide, z0, alpha);
      s0    = bsxfun(@rdivide, u,  alpha);
      

%[u,u,u,u]= fluid.pvt(initP, z0);

xrRef = initResSol (G, initP, s0, z0);
xrMs  = xrRef;
xwRef = initWellSol(W, initP);
xwMs  = xwRef;

gravity reset, gravity(useGrav);
porvol  = poreVolume(G, rock);
vd      = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);

solver_press = @(xr, xw, S, W, p0, dt) ...
   iterateMixedBO(xr, xw, G, rock, S,      ...
                           fluid, p0, dt, 'bc', bc, ...
                           'wells', W, 'src', src);
% solveMixedBO
% solveLinBO
% solveBlackOilWellSystem
% iterateMixedBO
%{
solver_pressMS = @(xr, xw, S, CS, W, p0, dt) ...
   solveCoarsePsysBO(xr, xw, G, CG, p, rock, S, CS, ...
                     fluid, p0, dt, 'bc', bc, 'wells', W, 'src', src);
%}
solver_pressMS = @(xr, xw, S, CS, W, p0, dt) ...
    iterateMixedBOMS(xr, xw, G, CG, p, rock, S, CS, fluid, ...
                 p0, dt, 'bc', bc, 'wells', W, 'src', src);
% solveLinBOMS                             
% iterateMixedBOMS

solver_sat   = @(xr, xw, tf, W) ...
   explicitBlackOil(G, xr, xw, tf, porvol, fluid, ...
                    'bc', bc, 'wells', W, 'src', src);

%--------------------------------------------------------------------------
%% Finish model initialisation --------------------------------------------

% Compute time (tau) until steady state and pressure time step.
%
[mu, mu, mu, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
tau = mu(1,end) .* rock.poro(1) * physDims(1)^2 / ...
      (rock.perm(1,1) * (pHi - pLo));
DT  = tau / 100;

% Construct basic pressure system components.
%
s   = bsxfun(@rdivide, u, u*ones([size(u,2), 1]));
mob = bsxfun(@rdivide, fluid.relperm(s), mu);
mob = mob * ones([size(mob,2), 1]);
S   = computeMimeticIP(G, rock, 'Verbose', verbose, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');
CS  = generateCoarseSystem(G, rock, S, CG, mob, 'bc', bc, 'src', src, 'BasisWeighting', 'poros');
W   = assembleWellSystem(G, W);
%WC  = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, 'src', src);

%--------------------------------------------------------------------------
%% Run model until convergence --------------------------------------------
%
T = 0;

pressRef      = zeros([401, G.cells.num]);
press         = zeros([401, G.cells.num]);
velRef        = zeros([401, numel(cellno)]);
vel           = zeros([401, numel(cellno)]);
pressRef(1,:) = xrRef.cellPressure;
press   (1,:) = xrRef.cellPressure;
velRef  (1,:) = xrRef.cellFlux;
vel     (1,:) = xrRef.cellFlux;

vdRef        = vd(DT, u);
vdMs         = vd(DT, u);

rhoS = fluid.surfaceDensity;

nonlin = @(xr, xw, solver) ...
   succSubst(xr, xw, solver, 'Tol', 1.0e-5/day, ...
             'MaxIt', 40, 'Verbose', true);

        
ep    = zeros(numSteps, 1);
emRef = zeros(numSteps, 1);
emMs  = zeros(numSteps, 1);
its   = zeros(numSteps, 1);

for step = 1 : numSteps
    disp(['------------------------- STEP: ' num2str(step) ' ---------------------'])
   S.RHS.volume_discrepancy = vdRef;
   % Compute reference solutions
   %
   solve_p = @(xr, xw) solver_press(xr, xw, S, W, xrRef.cellPressure, DT);
   [xrRef, xwRef] = nonlin(xrRef, xwRef, solve_p);

   %%{
   %overlap = 0;
   %xrMs0 = xrMs;
   %xwMs0 = xwMs;

   S.RHS.volume_discrepancy = vdMs;
   
   p0 = xrMs.cellPressure;
   solve_p = @(xr, xw) solver_pressMS(xr, xw, S, CS, WC, xrMs.cellPressure, DT);
   [xrMs, xwMs] = nonlin(xrMs, xwMs, solve_p);
   
   %useCorrection = false;
   xrMsUC = xrMs;
   for ii = 1:15
   if useCorrection
       [corrFunc] = solveMixedBOResiduals(xrMs, xwMs, G, rock, S, CG, ...
                                          fluid, p0, DT, ...
                                          'bc', bc, 'wells', W, 'src', src, ...
                                          'overlap', 20, 'tol', 1e-6);                                      
       for k = 1:numel(corrFunc)
           xrMs.cellPressure = xrMs.cellPressure + corrFunc(k).resSol.cellPressure;
           xrMs.cellFlux = xrMs.cellFlux + corrFunc(k).resSol.cellFlux;
           xrMs.facePressure = xrMs.facePressure + corrFunc(k).resSol.facePressure;
       end
       xrMs.faceFlux = cellFlux2faceFlux(G, xrMs.cellFlux);    
   end
   end
   
   %{
   itCount = 0;
   if useCorrection
       [xrMs, xwMs, corP, corV, itCount] = msCorrection(xrMs, xwMs, p0, G, S, ...
                                 W, rock, p, CG, CS, DT, fluid, 'bc', bc, ...
                                 'correctionOverlap', correctionOverlap, ...
                                 'tol', mscTol, 'mscIterNum', mscIterNum);
       xrMs.cellPressure = xrMs.cellPressure + sum(corP,2);
       xrMs.cellFlux = xrMs.cellFlux + sum(corV, 2);
       xrMs.faceFlux = cellFlux2faceFlux(G, xrMs.cellFlux);                                
   end
   its(step) = itCount;

   % --------------------------------------------------------------------- 
   %}   
FontSz = 14;

  
  subplot(3,2,[1 3]),
      plot(X, convertTo([xrRef.cellPressure, xrMs.cellPressure], barsa), ...
           'LineWidth', 2)
      v = axis; axis([v(1:2), -1, convertTo(pHi, barsa)])
      title('Cell Pressure', 'FontSize', FontSz), legend('Reference', 'MsMFEM'), set(gca, 'FontSize', FontSz);
   
   

   xrRef = solver_sat(xrRef, xwRef, DT, W );
   
   
   xrMs  = solver_sat(xrMs , xwMs , DT, WC);
   
   
   %xrMs = xrRef;
   %}

   [u, u, u, u] = fluid.pvt(xrRef.cellPressure, xrRef.z);
   vdRef        = vd(DT, u);

   [u, u, u, u] = fluid.pvt(xrMs.cellPressure, xrMs.z);
   vdMs         = vd(DT, u);

   massRef = bsxfun(@times, xrRef.z, rhoS);
   %%{
   massMs  = bsxfun(@times, xrMs .z, rhoS);
  
   relPDiff = (xrMs.cellPressure - xrRef.cellPressure)/(max(xrRef.cellPressure) - min(xrRef.cellPressure));
   ep(step) = norm(relPDiff, inf);
   %subplot(1,2,2)
   %plot(X, [massRef(:,1) massMs(:,1)], 'LineWidth', 2)

   FontSz = 14;
   subplot(3,2,5),
      %fp = convertTo(relPDiff, barsa);
      plot(X,relPDiff, 'LineWidth', 2)
      title(['Relative difference in cell pressure'])
      set(gca, 'FontSize', FontSz);
    

   
   subplot(3,2,6),
   relMassErrorRef = DT * (vdRef./porvol);
   relMassErrorMs  = DT * (vdMs./porvol);
   
   emRef(step)  = sum(abs(relMassErrorRef))/Nc;
   emMs(step)   = sum(abs(relMassErrorMs))/Nc;
   
   cla
   hold on
  plot(X, xrRef.s, 'LineWidth', 2); title('Saturations');
  plot(X, xrMs.s, '--', 'LineWidth', 2);
      %set(gca, 'XTick', []),
      set(gca, 'FontSize', FontSz);

   subplot(3,2,2),
      ff = convertTo([xrRef.faceFlux, xrMs.faceFlux], meter^3/day);
      plot(X,ff(1 : Nc, :)*[1 -1]'/norm(ff(1 : Nc, 1)), 'LineWidth', 2)
      title('Total flux - diff'), 
      %set(gca, 'XTick', []),
      set(gca, 'FontSize', FontSz);


   subplot(3,2,4),
   cla, hold on
      plot(X, [massRef], 'LineWidth', 2);
      plot(X, [massMs], '--','LineWidth', 2);
      title('Total mass');
      set(gca, 'FontSize', FontSz);
   %}
   T = T + DT;

   drawnow
    pressRef(step+1,:) = xrRef.cellPressure;
   velRef  (step+1,:) = xrRef.cellFlux;
   %%{
   press   (step+1,:) = xrMs .cellPressure;
   vel     (step+1,:) = xrMs .cellFlux;
end

