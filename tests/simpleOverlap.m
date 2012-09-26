%% Multiscale Pressure Solver with Overlap: 
% Use overlap in multiscale basis function to improve accuracy of the
% simulation.  We compare the fine-grid pressure solver and multiscale
% pressure solvers with and without overlap by solving the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a Cartesian grid with anisotropic, homogeneous permeability. 
%
% Overlap means that we extend the support of the multiscale basis
% functions, that is, we include more fine cells when computing the basis
% function. The use of overlap is particularly attractive in cases of
% highly heterogeneous reservoirs and wells that are placed near the corner
% of a coarse well-block. In the following, let the dashed lines be the
% coarse grid, the dotted lines represent the overlap, and let '*'
% represent a well.
%
%  Overlap in the reservoir basis function for a coarse face ij means that
% we include fine cells around the the neighboring coarse blocks when we
% compute the basis for face ij:
%               
%                  . . . . . . . . . .                           
%                . o-----------------o .
%                . |        |        | . 
%                . |   i    |   j    | .
%                . |        |        | .
%                . o-----------------o .  
%                  . . . . . . . . . .          
%                    
%  The flow near a well is described by a multiscale well basis. This
% basis normally only has support in the coarse well-block. If the well
% lies close to the corner of the well-block, this will be a poor
% description of the flow near the well because it will not describe the
% flow over the corner. In such cases we can use overlap around
% the well ('overlapWell'): 
%
%                           . . . . . 
%                           .       .
%                      o--------o   .
%                      |       *|   . 
%                      |        | . .
%                      |        |
%                      o--------o
% 
%  In some cases (e.g in case of highly hetereogeneous permeability) it
% might be necessary to use overlap around the entire well block
% ('OverlapBlock'):
%
%                    . . . . . . . . 
%                    . o---------o .
%                    . |         | . 
%                    . |    *    | .
%                    . |         | .
%                    . o---------o .
%                    . . . . . . . . 
%
%  This example shows the use of overlap around the well.   
%

%% Define and visualize the model
nx = 40; ny = 40; nz = 2;
cellDims = [nx, ny, nz];

G = cartGrid(cellDims);
G = computeGeometry(G);
diag_k     = [500, 500, 50] .* (1e-3*darcy());
rock.perm  = repmat(diag_k, [G.cells.num, 1]);
fluid = initSingleFluid();

% Set two vertical wells.
injCells  = [ 8,  8];
prodCells = [33, 33];

W = struct([]);
W = verticalWell(W, G, rock, injCells(1), injCells(2), 1:2, ...
                 'Type', 'rate', 'Val', 1*meter^3/day,      ...
                 'Radius', 0.1);
W = verticalWell(W, G, rock, prodCells(1), prodCells(2), 1:2, ...
                 'Type', 'bhp', 'Val', 0, 'Radius', 0.1);
W = assembleWellSystem(G, W);

% Visualize the model.
clf
plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [0.65, 0.65, 0.65]);
plotWell(G, W, 'radius', 0.1, 'color', 'r');
view(3), axis equal tight off



%% Partition the grid and assemble linear systems 
part  = partitionCartGrid(cellDims, [5, 5, 1]);
CG    = generateCoarseGrid(G, part, 'Verbose', true);

gravity off
S  = computeMimeticIP(G, rock, 'Verbose', true);

% Note: Only use overlap around wells, not in reservoir basisfunctions
CS = generateCoarseSystem (G, rock, S, CG, ones([G.cells.num,1]), ...
                           'Verbose', true, 'Overlap', 0);

%% Generate coarse well system with overlap around the well
% For overlap in the well basis we has two different choices: 'OverlapWell'
% and 'OverlapBlock'. The first includes cells around the well and the
% latter cells around the well-block (the coarse block containing the
% well). We apply overlap around the well for 'W1' and use no overlap for
% 'W' as a reference.                        
W1 = W;

mob = fluid.mob(struct('s', ones([G.cells.num, 1])));

% First no overlap
W  = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, ...
                              'OverlapWell',  0, 'OverlapBlock', 0);
% Then overlap around well
W1 = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W1, ...
                              'OverlapWell', 10, 'OverlapBlock', 0);

%% Solve the global flow problems 
%  The mass matrix B in the linear system will not be block-diagonal when
%  we use overlap. Thus, the lineary system system can not be reduced by
%  Schur-complement reduction as we do in the hybrid solver. Therefore, the
%  coarse system with overlap must be solved using option 'Solver', 'mixed'
%  as input to solveIncompFlowMS.

% Fine scale reference solution:
refR  = initResSol( G, 0.0);
refW  = initWellSol(W, 0.0);
[refR, refW] = solveIncompFlow(refR, refW, G, S, fluid, 'wells', W);

% Coarse scale - no overlap:
sR    = initResSol (G, 0.0);
sW    = initWellSol(W, 0.0);
[  sR,   sW] = solveIncompFlowMS(sR,  sW,  G, CG, part, S, CS, fluid, ...
                                 'wells', W, 'Solver', 'mixed');
% Coarse scale - overlap:
sR1   = initResSol( G, 0.0);
sW1   = initWellSol(W, 0.0);
[ sR1,  sW1] = solveIncompFlowMS(sR1, sW1, G, CG, part, S, CS, fluid, ...
                                 'wells', W1, 'Solver', 'mixed');

% Report pressure drop between wells.                         
dp = @(ws) num2str(convertTo(ws(1).pressure, barsa()));
disp(['DeltaP - Fine:          ', dp(refW)])
disp(['DeltaP - Ms no overlap: ', dp(sW  )])
disp(['DeltaP - Ms overlap:    ', dp(sW1 )])


%% Plot solution 
cellNo    = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
diverg    = @(v) accumarray(cellNo, abs(convertTo(v, meter^3/day)));
plot_flux = @(x) plotCellData(G, diverg(x));
well_bas  = @(w) S.BI * sparse(w.CS.basis{1}{1}, 1, w.CS.basis{1}{3}, ...
                               size(G.cellFaces,1), 1);
figure;
subplot(2,3,1)
   plot_flux(refR.cellFlux); title('Flux - fine scale'); cx = caxis;

subplot(2,3,2)
   plot_flux(sR.cellFlux)  ; title('Flux - no overlap'); caxis(cx)

subplot(2,3,3)
   plot_flux(sR1.cellFlux) ; title('Flux - overlap')   ; caxis(cx)

subplot(2,3,5)
   plot_flux(well_bas(W(1)));
   axis([0, 20, 0, 20]);  cx = caxis;
   title('Well basis - no overlap')

subplot(2,3,6)
   plot_flux(well_bas(W1(1)));
   axis([0, 20, 0, 20]); caxis(cx);
   title('Well basis - overlap')
