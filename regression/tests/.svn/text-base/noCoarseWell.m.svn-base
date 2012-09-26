function normDp = noCoarseWell(doPlot)
% compare relative pressure drop between wells for mimetic and multiscale
% where the coarse grid is equal to the fine grid. Compare pressure drop.
%
% RETURNS:
% normDp 

if nargin == 0
    doPlot = true;
end


% Geometry and rock type
cellDims = [10, 10, 10];
G = cartGrid(cellDims, cellDims);
G = computeGeometry(G);

rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

S = assembleMimeticSystem(G, rock);

% Wells 
W = struct([]);
W = verticalWell(W, G, rock, 10, 10, 1:10, 'Type', 'rate', ...
                 'Val', 1*meter^3/day, 'Radius', .1, 'Name', 'I');
W = addWell(G, rock, W, 1:10, 'Type','bhp', 'Val', 0, ...
            'Radius', .1, 'Dir', 'x', 'Name','P');
        
% Coarse and fine system:
p  = partitionUI(G, [10, 10, 10]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1])); 

W = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);


fluid = initSimpleFluid();
xrRef = initResSol(G, 0.0);
xwRef = initWellSol(W, 0.0);
xrMs  = initResSol(G, 0.0);
xwMs  = initWellSol(W, 0.0);

[xrRef, xwRef] = solveIncompFlow(  xrRef, xwRef, G, S, fluid, 'wells', W);
[xrMs , xwMs ] = solveIncompFlowMS(xrMs, xwMs, G, CG, p, S, CS, ...
                                   fluid, 'wells', W);
     
normDp=norm(1 - xwMs(1).pressure/xwRef(1).pressure);

% plot output
if doPlot
   subplot(2,2,1)
    plotGrid(G, [], 'FaceColor', 'none', 'EdgeColor', [.6, .6, .6]);
    plotGrid(G, W(1).cells, 'FaceColor', 'b', 'EdgeColor', 'b');
    plotGrid(G, W(2).cells, 'FaceColor', 'r', 'EdgeColor', 'r');
    title('Well cells')
    view(70, 15); axis tight; axis equal;

   subplot(2,2,2)
    plot(convertTo(xwRef(2).flux, meter^3/day), 'b'); hold on
    plot(convertTo(xwMs (2).flux, meter^3/day), 'r');
    legend('Fine','Multiscale')
    title('Producer inflow profile')
    
   subplot(2,2,3)
    plotCellData(G, convertTo(xrRef.cellPressure, barsa));
    title('Pressure Fine [bar]')
    view(70,15); axis tight; axis equal;
    cax = caxis; colorbar
    
   subplot(2,2,4)
    plotCellData(G, convertTo(xrMs.cellPressure, barsa));
    title('Pressure Coarse [bar]')
    view(70,15); axis tight; axis equal;
    caxis(cax); colorbar
end
