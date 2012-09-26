%% Test of use of global information in multiscale basis functions. 
% Compare the fine-grid and the multiscale pressure solver by solving the
% single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=-\frac{K}{\mu}\nabla p,$$
%
% for a Cartesian grid with isotropic, homogeneous permeability


nx = 20;
ny = 20;
nz = 1;

cellDims  = [nx, ny, nz];
verbose   = true;
G         = cartGrid(cellDims, cellDims);
G         = computeGeometry(G);
rock.perm = repmat(100*milli*darcy, [G.cells.num, 1]);
fluid     = initSimpleFluid();

% Set two wells, one vertical and one horizontal
W = struct([]);
W = verticalWell(W, G, rock, nx, ny, nz, ...
                  'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', .1, 'Name', 'I');
W = addWell(G, rock, W, 1, 'Type','bhp', ...
             'Val', 0, 'Radius', .1, 'Dir', 'x', 'Name', 'P');
          
src = [];
%src = addSource(src, [1, G.cells.num], [1 -1]./ day());
bc  = [];
%bc = addBc(bc, [1 11], 'flux', [0.5 -0.5]./day());


%% Set up solution structures
% Here we need four solution structures, two for each simulator to hold the
% solutions on the grid and in the wells, respectively.
xrRef = initResSol( G, 0.0);
xwRef = initWellSol(W, 0.0);
xrMs  = initResSol( G, 0.0);
xwMs  = initWellSol(W, 0.0);
xrMs2  = initResSol( G, 0.0);
xwMs2  = initWellSol(W, 0.0);

%% Partition the grid
% We partition the fine grid into a regular 5-by-5-by-2 coarse grid in
% index space so that each coarse block holds 8-by-8-by-5 fine cells. The
% resulting vector <p> has one entry per fine-grid cell giving the index of
% the corresponding coarse block. After the grid is partitioned in index
% space, we postprocess it to make sure that all blocks consist of a
% connected set of fine cells. This step is superfluous for Cartesian
% grids, but is required for grids that are only logically Cartesian (e.g.,
% corner-point and other mapped grids that may contain inactive or
% degenerate cells).
p  = partitionUI(G, [4, 4, 1]);
p  = processPartition  (G, p, 'Verbose', verbose);

CG = generateCoarseGrid(G, p, 'Verbose', verbose);

S  = computeMimeticIP(G, rock, 'Verbose', verbose);

W = assembleWellSystem(G, W);

[xrRef, xwRef] = solveIncompFlow(  xrRef, xwRef, G, S, ...
                                   fluid, 'wells', W, 'bc', bc, 'src', src,'Solver', 'hybrid');
                                
                                
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), 'global_inf', xrRef.faceFlux, ...
                           'Verbose', verbose,'src', src, 'bc', bc);
CS2 = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), ...
                           'Verbose', verbose, 'src', src,'bc', bc);

                        
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);
                                
                                
[xrMs , xwMs ] = solveIncompFlowMS(xrMs, xwMs, G, CG, p, S, CS, ...
                                   fluid, 'wells', W,'bc', bc, 'src', src, 'Solver', 'hybrid');

                                
[xrMs2 , xwMs2 ] = solveIncompFlowMS(xrMs2, xwMs2, G, CG, p, S, CS2, ...
                                   fluid, 'wells', W,'bc', bc, 'src', src, 'Solver', 'hybrid');
                                
                                
dp = @(x) num2str(convertTo(x(1).pressure, barsa));
disp(['DeltaP - Fine:     ', dp(xwRef)]);
disp(['DeltaP - Ms_global:', dp(xwMs )]);
disp(['DeltaP - Ms_old:   ', dp(xwMs2 )]);

disp(['rel_norm_flux(ref-Ms_global): ',num2str(norm(xrRef.faceFlux-xrMs.faceFlux)/norm(xrRef.faceFlux))])
disp(['rel_norm_flux(ref-Ms_old):    ',num2str(norm(xrRef.faceFlux-xrMs2.faceFlux)/norm(xrRef.faceFlux))])
disp(['rel_norm_pres(ref-Ms_global): ',num2str(norm(xrRef.cellPressure-xrMs.cellPressure)/norm(xrRef.cellPressure))])
disp(['rel_norm_pres(ref-Ms_old):    ',num2str(norm(xrRef.cellPressure-xrMs2.cellPressure)/norm(xrRef.cellPressure))])
cellNo    = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
plot_var  = @(x) plotCellData(G, x);
plot_flux = @(x) plot_var(accumarray(cellNo, abs(convertTo(x.cellFlux, meter^3/day))));

%% plot output
f = figure;
subplot(3,2,1)
   plot_flux(xrRef);
   view(3), camproj perspective, axis tight equal, camlight headlight
   title('Reference solution')
   cax2 = caxis; colorbar
subplot(3,2,3)
   plot_flux(xrMs);
   caxis(cax2); colorbar
   title('MS with global info')
   view(3), camproj perspective, axis tight equal, camlight headlight
subplot(3,2,5)
   plot_flux(xrMs2);
   caxis(cax2); colorbar
   title('MS without global info')
   view(3), camproj perspective, axis tight equal, camlight headlight
subplot(3,2,2)
   plotCellData(G, convertTo(xrRef.cellPressure, barsa));
   title('Pressure Fine')   
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax = caxis; colorbar

subplot(3,2,4)
   plotCellData(G, convertTo(xrMs.cellPressure, barsa));
   title('Pressure Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax); colorbar
subplot(3,2,6)
   plotCellData(G, convertTo(xrMs2.cellPressure, barsa));
   title('Pressure Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax); colorbar   
   
   
