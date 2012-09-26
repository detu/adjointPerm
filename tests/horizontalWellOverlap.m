% horizonalWellOverlap;
% Long horizontal well using block- and well overlap

% Clear
cellDims = [60, 60, 5];
resDims  = [500 500 20];
verbose  = true;
sizLayer = prod(cellDims(1:2));

blockOverlap = 6;
wellOverlap  = 6;

G = computeGeometry(cartGrid(cellDims, resDims));

% layered perm-field
lo = ones([sizLayer, 1]);
lp = exp(5*( rand(cellDims(3),1) - .5) );
rock.perm = kron(lp, lo) * milli*darcy;

% 3 vert inj, 1 hor prod
W = addWell(G, rock, [], 49*60+55 : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf');
W = addWell(G, rock, W, 1         : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf');
W = addWell(G, rock, W, 30        : sizLayer : G.cells.num,      ...
            'Type', 'rate', 'Val', 1/3*meter^3/day, 'Radius', 1, ...
            'InnerProduct', 'ip_tpf');

W = addWell(G, rock, W, 2*sizLayer + 14*60 + (10:55),         ...
            'Type', 'bhp', 'Val', 0, 'Radius', 1, 'Dir', 'x', ...
            'InnerProduct', 'ip_tpf');

part  = partitionUI(G, [4, 4, 1]);
part2 = processPartition  (G, part,  'Verbose', verbose);
CG    = generateCoarseGrid(G, part2, 'Verbose', verbose);

S = computeMimeticIP(G, rock,            ...
                           'Verbose', verbose, ...
                           'Type', 'comp_hybrid',   ...
                           'InnerProduct', 'ip_tpf');

CS = generateCoarseSystem (G, rock, S, CG, ones([G.cells.num, 1]), ...
                           'Verbose', verbose, 'Overlap', blockOverlap);


W = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), ...
                             rock, W, 'OverlapWell', wellOverlap,  ...
                             'OverlapBlock', blockOverlap);

xrRef = initResSol (G, 0.0);
xwRef = initWellSol(W, 0.0);
xrMs  = initResSol (G, 0.0);
xwMs  = initWellSol(W, 0.0);
fluid = initSimpleFluid();

[xrRef, xwRef] = solveIncompFlow  (xrRef, xwRef, G, S, fluid, ...
                                   'wells', W, 'Solver', 'mixed');
[xrMs,  xwMs ] = solveIncompFlowMS(xrMs, xwMs, G, CG, part2, S, CS, ...
                                   fluid, 'wells', W, 'Solver', 'mixed');
    
disp(['DeltaP - Fine: ', num2str(convertTo(xwRef(1).pressure, barsa))])
disp(['DeltaP - Ms:   ', num2str(convertTo(xwMs (1).pressure, barsa))])

%% plot output
f = figure;
subplot(2,2,1)
   plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.8, .8, .8])
   plotGrid(G, reshape([W(1:3).cells], [], 1), 'FaceColor', 'b', ...
            'EdgeColor', 'black');
   plotGrid(G, W(4).cells, 'FaceColor', 'r', 'EdgeColor', 'black');
   axis tight
   title('Well-cells')
   
subplot(2,2,2)
   plot(-convertTo(xwRef(4).flux, meter^3/day), 'b'); hold on
   plot(-convertTo(xwMs (4).flux, meter^3/day), 'r');
   axis([0 46 0 .1]); 
   legend('Fine', 'Multiscale')
   title('Producer inflow profile')

subplot(2,2,3)
   plotCellData(G, convertTo(xrRef.cellPressure, barsa));
   title('Pressure Fine [bar]')
   axis tight , cax = caxis; colorbar

subplot(2,2,4)
   plotCellData(G, convertTo(xrMs.cellPressure, barsa));
   title('Pressure Coarse [bar]')
   axis tight, caxis(cax); colorbar
