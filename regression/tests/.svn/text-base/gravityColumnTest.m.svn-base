function normP = gravityColumnTest(doPlot)
% Gravity example for multiscale: Put pressure bc on top of column and let
% pressure be governed by gravity. Coarse grid is equal to fine grid in
% z-direction. 
%
% RETURNS: 
% normP: norm of difference in cellPressure.

if nargin == 0
    doPlot = true;
end

gravity on;

G = cartGrid([2, 2, 10], [1, 1, 10]);
G = computeGeometry(G);

rock.perm  = repmat(100*darcy()/1000, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);

% Impose pressure bc
S  = assembleMimeticSystem(G, rock);
bc = pside([], G, 'TOP', 1:2, 1:2, 100e5); %, 'sat', repmat([1,0,0], [4,1]));

                     
% Coarse system:
p  = partitionUI(G, [2, 1, 10]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ones([G.cells.num, 1]), 'bc', bc); 

% Solve fine and coarse system 
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);

[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);
                                
                                
%Check multiscale pressure against fine-scale                                
normP = norm(xrMs.cellPressure - xrRef.cellPressure)/ ... 
             norm(xrRef.cellPressure);                        
                                                                       
% plot output
if doPlot
   figure;
    subplot(1, 2, 1)
    plotFaces(G, 1:G.faces.num, xrRef.facePressure);
    plotCellData(G, xrRef.cellPressure);
    title('Pressure')
    cax=caxis;
    view(3)
    colorbar
    
    subplot(1, 2, 2)
    % NB: uses cellPressure instead of facePressure
    plotCellData(G, xrMs.cellPressure);
    title('Pressure Coarse')
    view(3);
    caxis(cax)
    colorbar
end