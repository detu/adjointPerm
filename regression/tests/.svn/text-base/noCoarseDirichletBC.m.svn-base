function normFP = noCoarseDirichletBC(doPlot)
% Simulate system with dirichletBC.
% Use no coarsening, i.e grid for multiscale = grid for mimetic
% 
% RETRUNS:
% normFP = normFlux + normPressure(global), where
% normFlux = norm(sol_ms.faceFlux(activeBC) - sol_mim.faceFlux(activeBC))/ ...
%              norm(sol_mim.faceFlux(activeBC));
% normPressure = norm(sol_ms.cellPressure-sol_mim.cellPressure)/...
%                  norm(sol_mim.cellPressure)    

if nargin==0
    doPlot = true;
end

gravity off;

% Define geometry and rock type
cellDims = [10, 10, 1];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat([1000, 100, 10], [G.cells.num, 1]);
rock.poro = repmat(0.3,             [G.cells.num, 1]);

% Impose Dirichlet BCs.
facesBC1 = find(G.faces.tag == 1);
facesBC2 = find(G.faces.tag == 2);

S  = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'pressure', 300, 'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'pressure', 250, 'sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [10, 10, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]),...
    'bc',bc); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);
                          
[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);
                                             
% Calculate norm of flux on active boundary 
normFlux     = norm(xrMs.faceFlux([facesBC1; facesBC2]) - ...
                     xrRef.faceFlux([facesBC1; facesBC2]))/  ...
                     norm(xrRef.faceFlux([facesBC1; facesBC2]));
normFP       = normFlux;

% %% plot output
% f = figure;
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
    caxis(cax2); colorbar
end