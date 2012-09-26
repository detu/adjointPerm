function normDp = wellHetero2D(doPlot)
% simple example of cartesian 2D grid with permeability from SPE10 with one
% pressure and one rate controlled well. 
%
% RETURNS:
% normDp: norm of difference in pressure drop between wells. 
%
% COMMENTS:
% Should use porosity from SPE10

if nargin == 0
    doPlot = true;
end

% Define geometry and rock type
nx = 20; ny = 40; nz = 1;
cellDims = [nx, ny, nz];
G = cartGrid(cellDims);
G = computeGeometry(G);

load Kx;
startX = 20; startY = 60; startZ = 40;
rock.perm  = Kx(startX:startX+nx-1,startY:startY+ny-1,startZ:startZ+nz-1);
rock.perm = rock.perm(:);
rock.poro = repmat(0.3, [G.cells.num, 1]);
clear kx;

% Put wells
S = assembleMimeticSystem(G, rock, 'Gravity', false);
W = struct([]);
W = verticalWell(W, G, rock, nx, ny, 1:nz, 'Type', 'rate', ...
                 'Val', 1*meter^3/day, ...
                 'Radius', .1, 'Name',  'I');

W = verticalWell(W, G, rock, 1, 1, 1:nz, 'Type', 'bhp', ...
                 'Val', 0, 'Radius', .1, 'Name',  'P');

% Coarse system:
p  = partitionUI(G, [5, 10, 1]);
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

[xrRef, xwRef] = solveIncompFlow  (xrRef, xwRef, G, S, fluid, 'wells', W);
[xrMs,  xwMs ] = solveIncompFlowMS(xrMs, xwMs, G, CG, p, S, CS, ...
                                   fluid, 'wells', W);     

normDp = norm(1 - xwMs(1).pressure/xwRef(1).pressure);

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
