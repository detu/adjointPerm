function normFP = dirichletNeumannBC(doPlot)
% Impose dirichlet bc and neumann bc on left and right side of a simple
% system. Compare fine scale flux and multiscale flux over the pressure
% boundary faces and mean fince scale pressure against mean coarse pressure
% over flux boundary faces. 
%
% RETURNS: 
%
% normFP = normFlux + normPressure
%

if nargin == 0
    doPlot = true;
end

% Define geometry and rock type
cellDims = [50, 50, 4];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat([1000, 100, 10], [G.cells.num, 1]);
rock.poro = repmat(0.3,             [G.cells.num, 1]);

% Impose Dirichlet BCs and put empty well.
facesBC1 = find(G.faces.tag==1);
facesBC2 = find(G.faces.tag==2);

S = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'pressure',  1, 'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'flux',      1, 'sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [5, 5, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num,1]), 'bc', bc); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);

[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);

% Calculate norm of flux on active boundary 
normFlux     = norm(xrMs.faceFlux(facesBC1) - xrRef.faceFlux(facesBC1))/ ...
                     norm(xrRef.faceFlux(facesBC1));
normFP       = normFlux;
                       
% plot output
if doPlot
    subplot(2,2,1)
    plotCellData(G, xrRef.cellPressure);
    title('Pressure Fine')
    view(3); axis tight; axis equal;
    cax = caxis; colorbar
    
    subplot(2,2,2)
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    view(3); axis tight; axis equal;
    caxis(cax);
    colorbar
    
    subplot(2,2,3)
    title('Flux intensity Fine')
    plotCellData(G, sqrt(S.C .' * abs(xrRef.cellFlux)));
    view(3); axis tight; axis equal;
    cax2 = caxis; colorbar
    
    subplot(2,2,4)
    title('Flux intensity Coarse')
    plotCellData(G, sqrt(S.C .' * abs(xrMs.cellFlux)));
    view(3); axis tight; axis equal;
    caxis(cax2);
    colorbar
end
    
