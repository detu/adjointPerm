% simpleMSBCExample
% The example shows how to put Dirichlet boundary conditions on the system.

% $Id: simpleMSBCExample.m 1916 2009-03-27 09:32:44Z jrn $

clear
cellDims = [40, 40, 10];
verbose  = true;

gravity off;

G = cartGrid(cellDims, cellDims);
G = computeGeometry(G);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3            , [G.cells.num, 1]);
W = struct([]);
W = verticalWell(W, G, rock, 40, 40, 1:10, ...
                 'Type', 'rate', 'Val', 1*meter^3/day, ...
                 'Radius', .1, 'Name', 'I');
W = addWell(G, rock, W, 1:40, 'Type','bhp', ...
            'Val', 0, 'Radius', .1, 'Dir', 'x', 'Name', 'P');

fluid = initSimpleFluid();
xrRef = initResSol(G, 0.0);

p  = partitionUI(G, [5, 5, 2]);
p  = processPartition  (G, p, 'Verbose', verbose);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);
S  = computeMimeticIP(G, rock, 'Verbose', verbose);
%
bc = pside([], G, 'LEFT', 1:cellDims(2), 1:cellDims(3), 0);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), ...
                          'Verbose', verbose, 'bc', bc);
     
W = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);

[xrRef, xwRef] = solveIncompFlow  (xrRef, initWellSol(W, 0.0), G, S, fluid, ...
                                   'bc', bc, 'wells', W);
[xrMs , xwMs ] = solveIncompFlowMS(initResSol(G, 0.0), initWellSol(W, 0.0), ...
                                   G, CG, p, S, CS, fluid, 'wells', W, ...
                                   'bc', bc);

%% plot output
f = figure;
cellNo = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
plot_var = @(x) plotCellData(G, x);
plot_pres = @(x) plot_var(convertTo(x.cellPressure, barsa));
plot_flux = @(x) plot_var(log10(accumarray(cellNo, ...
                          abs(convertTo(x.cellFlux, meter^3/day)))));
subplot(2,2,1)
   plot_pres(xrRef); title('Pressure Fine [bar]')
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax = caxis; colorbar

subplot(2,2,2)
   plot_pres(xrMs); title('Pressure Coarse [bar]')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax); 
   colorbar
   
subplot(2,2,3)
   plot_flux(xrRef); title('Flux intensity Fine')
   view(3), camproj perspective, axis tight equal, camlight headlight
   cax2 = caxis; colorbar

subplot(2,2,4)
   plot_flux(xrMs); title('Flux intensity Coarse')
   view(3), camproj perspective, axis tight equal, camlight headlight
   caxis(cax2); colorbar
