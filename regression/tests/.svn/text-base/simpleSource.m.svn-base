function normFlux = simpleSource(doPlot)
% Add two sources to system and look at norm of difference in flux 
%
% RETURNS: 
% normFlux

if nargin == 0
    doPlot = true;
end

gravity off;

% Define geometry and rock type
cellDims = [40, 40, 1];
G = cartGrid(cellDims);
G = computeGeometry(G);

%rock.perm  = repmat([1000, 100, 10], [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Put sources.
S = assembleMimeticSystem(G, rock);
src = addSource([], [4*40+4 G.cells.num-(4*40+4)], [-1 1], ...
                'sat', repmat([1,0,0], [2,1]));

% Coarse system:
% Hack to support partitioning function 'partitionUI'.
[ijk{1:3}] = ind2sub(cellDims, (1 : G.cells.num).');
G.cells.ijkMap = [ijk{:}];

p  = partitionUI(G, [5, 5, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), 'src', src); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);

[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'src', src);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'src', src);

% calculate norm of difference in pressure drop                             
normFlux = norm(xrMs.faceFlux - xrRef.faceFlux)/norm(xrRef.faceFlux);  

%% plot output
if doPlot
    subplot(2,2,1)
    plotCellData(G, xrRef.cellPressure);
    title('Pressure Fine')
    axis tight; axis equal;
    cax = caxis; colorbar
    
    subplot(2,2,2)
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    axis tight; axis equal;
    caxis(cax);
    colorbar
    
    subplot(2,2,3)
    title('Flux intensity Fine')
    plotCellData(G, sqrt(S.C .' * abs(xrRef.cellFlux)));
    axis tight; axis equal;
    cax2 = caxis; colorbar
    subplot(2,2,4)
    title('Flux intensity Coarse')
    plotCellData(G, sqrt(S.C .' * abs(xrMs.cellFlux)));
    axis tight; axis equal;
    caxis(cax2);
    colorbar
end