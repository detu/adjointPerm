function normDp = neumannBC(doPlot)
% Impose neumann bc on simple system. Compare fine scale pressure drop and
% multiscale pressure drop by looking at the relative error norm the
% mean pressure over the active boundary faces.  
%
% RETURNS: 
% normDp = norm(1-Dp_ms/Dp_ref))    

if nargin == 0
    doPlot = true;
end

gravity off;

% Define geometry and rock type
cellDims = [50, 10, 2];
G = cartGrid(cellDims);
G = computeGeometry(G);

%rock.perm  = repmat([1000, 100, 10], [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Impose Neumann BCs.
% active boundary faces:
facesBC1 = find(G.faces.tag == 1);
facesBC2 = find(G.faces.tag == 2);

S  = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'flux',  1, 'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'flux', -1, 'sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [5, 1, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num,1]), 'bc', bc); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);
                          
[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);                                

% calculate norm of difference in pressure drop
% find active coarse boundary faces
cFacesBC1 = (CG.faces.tag == 1);
cFacesBC2 = (CG.faces.tag == 2); 

normFlux = norm(xrMs.faceFlux([facesBC1; facesBC2]) - ... 
                    xrRef.faceFlux([facesBC1; facesBC2]))/...
                    norm( xrRef.faceFlux([facesBC1; facesBC2]));                                           

%{
Dp_fine = abs(mean(xrRef.facePressure(facesBC1)) - ...
              mean(xrRef.facePressure(facesBC2)));
Dp_ms   = abs(mean(xrMs.blockFacePressure(cFacesBC1)) - ...
              mean(xrMs.blockFacePressure(cFacesBC2)));
normDp  = norm(1-Dp_ms/Dp_fine);                                 
%}
normDp  = normFlux;

% plot output
if doPlot
    subplot(2,2,1)
    plotCellData(G, xrRef.cellPressure);
    title('Pressure Fine')
    view(3); camproj perspective; axis tight; axis equal; camlight headlight;
    cax = caxis; colorbar
    
    subplot(2,2,2)
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    view(3); camproj perspective; axis tight; axis equal; camlight headlight;
    caxis(cax);
    colorbar
    
    subplot(2,2,3)
    title('Flux intensity Fine')
    plotCellData(G, sqrt(S.C .' * abs(xrRef.cellFlux)));
    
    view(3);camproj perspective; axis tight; axis equal; %camlight headlight;
    cax2 = caxis; colorbar
    subplot(2,2,4)
    title('Flux intensity Coarse')
    
    plotCellData(G, sqrt(S.C .' * abs(xrMs.cellFlux)));
    view(3);camproj perspective; axis tight; axis equal; %camlight headlight;
    caxis(cax2);
    colorbar
end
