function normDp = simpleWell(doPlot)
% simple example of cartesian grid with unit permeability with one pressure
% and one rate controlled well. 
%
% RETURNS:
% normDp: norm of difference in pressure drop between wells. 
%

if nargin == 0
    doPlot = true;
end

% Define geometry and rock type
cellDims = [40, 40, 10];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3,             [G.cells.num, 1]);

% Put wells
S = assembleMimeticSystem(G, rock);
W = struct([]);
W = verticalWell(W, G, rock, 40, 40, 1:10, 'Type', 'rate', ...
                'Val', 1*meter^3/day, 'Radius', .1, 'Name', 'I');
W = addWell(G, rock, W, 1:40, 'Type','bhp', 'Val', 0, 'Radius', .1, ...
            'Dir', 'x', 'Name','P');
         
% Coarse system:
p  = partitionUI(G, [5, 5, 2]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1])); 

W = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);

% Solve fine and coarse system 
fluid = initSimpleFluid();
xrRef = initResSol(G, 0.0);
xwRef = initWellSol(W, 0.0);
xrMs  = initResSol(G, 0.0);
xwMs  = initWellSol(W, 0.0);

[xrRef, xwRef] = solveIncompFlow(  xrRef, xwRef, G, S, fluid, 'wells', W);
[xrMs , xwMs ] = solveIncompFlowMS(xrMs, xwMs, G, CG, p, S, CS, ...
                                   fluid, 'wells', W); 

normDp = norm(1 - xwMs(1).pressure/xwRef(1).pressure);

% plot output
if doPlot,
   cellNo     = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
   plot_var   = @(x) plotCellData(G, x);
   plot_press = @(x) plot_var(convertTo(x.cellPressure, barsa));
   plot_flux  = @(x) plot_var(accumarray(cellNo, ...
                              abs(convertTo(x.cellFlux, meter^3/day))));

   subplot(2,2,1)
   plot_press(xrRef); title('Pressure Fine [bar]')
   axis tight equal, view(3), cax = caxis; colorbar

   subplot(2,2,2)
   plot_press(xrMs); title('Pressure Coarse [bar]')
   axis tight equal, view(3), caxis(cax); colorbar

   subplot(2,2,3)
   plot_flux(xrRef); title('Flux intensity Fine')
   axis tight equal, view(3), cax2 = caxis; colorbar

   subplot(2,2,4)
   plot_flux(xrMs); title('Flux intensity Coarse')
   axis tight equal, view(3), caxis(cax2); colorbar
end
