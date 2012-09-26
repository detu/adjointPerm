function normFlux = periodicBC(doPlot)
% Impose a combination of periodic (on left (1) and right (2) side ) and
% no-flow (up (4) and down (3)) flux bc. Compare fine scale flux and
% multiscale flux over the entire simulation grid by looking at the norm.
%
% RETURNS: 
% normFlux 

if nargin == 0
    doPlot = true;
end

gravity off;
% Define geometry and rock type
cellDims = [20, 20, 1];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat(100, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Impose Dirichlet BCs.
facesBC1 = find(G.faces.tag==1);
facesBC2 = find(G.faces.tag==2);

S = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'flux', -(1:length(facesBC1))',   'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'flux',  (length(facesBC2):-1:1)','sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [5, 5, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), 'bc', bc); 

% Solve fine and coarse system fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);
fluid    = initSimpleFluid();
                                 
[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);

% Calculate norm of flux on active boundary 
normFlux = norm(xrMs.faceFlux - xrRef.faceFlux)/norm(xrRef.faceFlux);                                          

% plot output
if doPlot
    subplot(2,2,1)
    plotCellData(G, xrRef.cellPressure);
    title('Pressure Fine')
    axis equal; axis tight;
    cax = caxis; colorbar
    
    subplot(2,2,2)
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    axis equal; axis tight;
    caxis(cax); colorbar
    
    subplot(2,2,3)
    title('Flux intensity Fine')
    plotCellData(G, sqrt(S.C .' * abs(xrRef.cellFlux)));
    axis equal; axis tight;
    cax2 = caxis; colorbar
    
    subplot(2,2,4)
    title('Flux intensity Coarse')
    plotCellData(G, sqrt(S.C .' * abs(xrMs.cellFlux)));
    axis equal; axis tight;
    caxis(cax2);
    colorbar
end
