%g =  computeGeometry(cartGrid([20,10,10]));
g  = computeGeometry(cartGrid([6,6]));
g  = nodeFaces(g);
clear rock
for i = 1:size(g.nodes.coords, 2),
   rock.perm(:,i) = logNormLayers([g.cartDims, 1], 1);
end

q   = zeros(g.cells.num, 1); q(1) =1; q(end)=-1;
src = [];%addSource([], [1, g.cells.num], [1, -1]);
bc  = fluxside([], g, 'left',  1:g.cartDims(2), 1,  1);
bc  = fluxside(bc, g, 'right', 1:g.cartDims(2), 1, -1);

bc  = pside([], g, 'left',  1:g.cartDims(2), 1,  2);
bc  = pside(bc, g, 'right', 1:g.cartDims(2), 1,  1);

T1  = computeMultiPointTrans(g, rock);
xr1 = incompMPFA(initResSol(g, 0, 0), [], g, T1, initSingleFluid(),...
      'src', src, 'bc', bc);

T2  = computeTrans(g, rock);
xr2 = incompTPFA(initResSol(g, 0, 0), [], g, T2, initSingleFluid(), ...
      'src', src, 'bc', bc);

cellno = rldecode(1:g.cells.num, double(g.cells.numFaces), 2) .';
cf     = g.cellFaces(:,1);

clf,

subplot(2,2,1);
cf1 = accumarray(cellno, xr1.faceFlux(cf));
plotCellData (g, cf1);
colorbar
cax = caxis;
axis equal tight

subplot(2,2,2);
cf2 = accumarray(cellno, xr2.faceFlux(cf));
plotCellData (g, cf2);
colorbar
caxis(cax);
axis equal tight

subplot(2,2,3);
plotCellData (g, xr1.cellPressure);
colorbar
axis equal tight

subplot(2,2,4);
plotCellData (g, xr2.cellPressure);
colorbar
axis equal tight


fprintf('Difference in fluxes: \t\t%e\nDifference in cell pressures:\t%e\n', ...
    norm(cf1-cf2, inf), norm(xr1.cellPressure-xr2.cellPressure, inf));
