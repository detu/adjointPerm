function normFlux = dirichletBCHetero2D(doPlot)
% Impose dirichlet bc on system with cartesian 2D grid and permeability
% from SPE10 system. Compare fine scale flux and multiscale flux over the
% active boundary faces by looking at the norm.
%
% RETURNS:  
% normFlux  
%
% COMMENTS:
% Should use porosity from SPE10.

if nargin == 0
    doPlot = true;
end

gravity off;

% Define geometry and rock type
nx = 50; ny = 50; nz = 1;
cellDims = [nx, ny, nz];
G = cartGrid(cellDims);
G = computeGeometry(G);

load Kx;
startX = 10; startY = 10; startZ = 46; 
rock.perm  = Kx(startX:startX+nx-1,startY:startY+ny-1,startZ:startZ+nz-1);
rock.perm = rock.perm(:);
clear Kx;
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Impose Dirichlet BCs.
facesBC1 = find(G.faces.tag==3);
facesBC2 = find(G.faces.tag==4);

S = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'pressure', 300, 'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'pressure', 250, 'sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [5, 5, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), 'bc', bc); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);

[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);                                
                                
% Calculate norm of flux on active boundary 
normFlux = norm(xrMs.faceFlux([facesBC1; facesBC2]) - ... 
                    xrRef.faceFlux([facesBC1; facesBC2]))/...
                    norm( xrRef.faceFlux([facesBC1; facesBC2]));                                           

% plot output
if doPlot
    subplot(3,2,1)
    plotCellData(G, xrRef.cellPressure);
    title('Pressure Fine')
    axis equal;  axis tight; view(3)
    cax = caxis; colorbar
    
    subplot(3,2,2)
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    axis equal; axis tight; view(3);
    caxis(cax);
    colorbar
    
    subplot(3,2,3)
    title('Flux intensity Fine')
    plotCellData(G, sqrt(S.C .' * abs(xrRef.cellFlux)));
    axis equal; axis tight; view(3);
    cax2 = caxis; colorbar
    
    subplot(3,2,4)
    title('Flux intensity Coarse')
    plotCellData(G, sqrt(S.C .' * abs(xrMs.cellFlux)));
    axis equal; axis tight; view(3)
    caxis(cax2); colorbar
    
    subplot(3,2,5)
    title('Perm')
    plotCellData(G, log(rock.perm));
    axis equal; axis tight; view(3)
    outlineCoarseGrid(G, p);
    colorbar
end
