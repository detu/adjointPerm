% readSimGridExample


caseName = 'data/model3';
G = readSimData(caseName, [30 50 5]);

% --- reading perm should be done in above function ---
loc = [caseName, '-perm.dat'];
[fid, msg] = fopen(loc, 'r');
rock.perm = fscanf(fid, '%f', inf);

% -----------------------------------------------------


%% Assemble system
S = computeMimeticIP(G, rock, 'Verbose', true);

%% Introduce wells
wellCells1 = find( and(G.cells.ijkMap(:,1) == 5, G.cells.ijkMap(:,2) == 20) );
wellCells2 = find( and(G.cells.ijkMap(:,1) == 43, G.cells.ijkMap(:,2) == 10) );

radius = .1;

fluid   = initSimpleFluid();
resSol  = initResSol(G, 0.0);
wellSol = initWellSol(W, 0.0);

W = addWell(G, rock, [], wellCells1, 'Type','rate','Val', 1,'Radius', radius);
W = addWell(G, rock, W, wellCells2, 'Type','bhp' ,'Val', 0, 'Radius',radius);

W = assembleWellSystem(G, W);

%% Solve system
[resSol, wellSol] = solveIncompFlow(resSol, wellSol, G, S, fluid, ...
                                    'wells', W);

%% plot output
f = figure;
subplot(2,2,1)
plotGrid(G, [], 'FaceColor', 'none', 'EdgeColor', [.6 .6 .6]);
plotGrid(G, W(1).cells, 'FaceColor', 'b', 'EdgeColor', 'b');
plotGrid(G, W(2).cells, 'FaceColor', 'r', 'EdgeColor', 'r');
title('Well-cells')
view(3); camproj perspective; axis tight; axis off; camlight headlight;

subplot(2,2,2)
plotCellData(G, resSol.cellPressure); 
title('Pressure')
view(3); camproj perspective; axis tight; axis off; camlight headlight;

subplot(2,2,3)
plotCellData(G, sqrt(S.C'*abs(resSol.cellFlux)) );
title('Sqrt - Flux intensity')
view(3); camproj perspective; axis tight; axis off; camlight headlight;

subplot(2,2,4)
plot(-wellSol(2).flux)
title('Producer inflow profile')
