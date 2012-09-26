function normFlux = simpleWellStrip(doPlot)
% Add two wells to system and look at norm of difference in flux 
%
% RETURNS: 
% normFlux: norm of difference in faceFlux       

if nargin == 0
    doPlot = true;
end

% Define geometry and rock type
cellDims = [40, 21, 1];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Put wells
S = assembleMimeticSystem(G, rock);
W = struct([]);

W = addWell(G, rock, W, 10*40+4, 'Type','bhp', 'Val', 0, ...
            'Radius', .1, 'Dir', 'x', 'Name','P');
W = addWell(G, rock, W, G.cells.num-(10*40+4), 'Type', 'rate', ...
            'Val', 1*meter^3/day, 'Radius', .1, 'Dir', 'x', 'Name', 'P');
         
% Coarse system:
p  = partitionUI(G, [5, 3, 1]);
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

[xrRef] = solveIncompFlow(  xrRef, xwRef, G, S, fluid, 'wells', W);
[xrMs]  = solveIncompFlowMS(xrMs, xwMs, G, CG, p, S, CS, ...
                                   fluid, 'wells', W); 
                            
                                
% calculate norm of difference in flux                        
normFlux = norm(xrMs.faceFlux - xrRef.faceFlux)/norm(xrRef.faceFlux);

% plot output
if doPlot
   cellNo     = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
   plot_var   = @(x) plotCellData(G, x);
   plot_press = @(x) plot_var(convertTo(x.cellPressure, barsa));
   plot_flux  = @(x) plot_var(accumarray(cellNo, ...
                              abs(convertTo(x.cellFlux, meter^3/day))));

   subplot(2,2,1)
   plot_press(xrRef); title('Pressure Fine [bar]')
   axis tight, cax = caxis; colorbar

   subplot(2,2,2)
   plot_press(xrMs); title('Pressure Coarse [bar]')
   axis tight, caxis(cax); colorbar

   subplot(2,2,3)
   plot_flux(xrRef); title('Flux intensity Fine')
   axis tight, cax2 = caxis; colorbar

   subplot(2,2,4)
   plot_flux(xrMs); title('Flux intensity Coarse')
   axis tight, caxis(cax2); colorbar
end
