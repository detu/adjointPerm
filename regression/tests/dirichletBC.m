function normFlux = dirichletBC(doPlot)
% Impose dirichlet BC on simple system (put pressure on left and right
% side).  Compare fine scale flux and multiscale flux over the active
% boundary faces by looking at the relative error norm.
%
% RETURNS:
%    normFlux

if nargin == 0,
   doPlot = true;
end

% Define geometry and rock type
nx = 50;
ny = 50;
nz = 4;
cellDims = [nx, ny, nz];
G = cartGrid(cellDims);
G = computeGeometry(G);

rock.perm  = repmat([1000, 100, 10], [G.cells.num, 1]);
rock.poro = repmat(0.3,             [G.cells.num, 1]);

% Impose Dirichlet BCs.
facesBC1 = find(G.faces.tag==1);
facesBC2 = find(G.faces.tag==2);

S  = assembleMimeticSystem(G, rock);
bc = addBC([], facesBC1, 'pressure', 300, 'sat', [1,0,0]);
bc = addBC(bc, facesBC2, 'pressure', 250, 'sat', [1,0,0]);

% Coarse system:
p  = partitionUI(G, [5, 5, 1]);
p  = processPartition  (G, p);
CG = generateCoarseGrid(G, p);
CS = generateCoarseSystem(G, rock, S, CG, ...
                          ones([G.cells.num,1]), 'bc', bc);

% Solve fine and coarse system
fluid    = initSimpleFluid();
xrRef    = initResSol(G, 0.0);
xrMs     = initResSol(G, 0.0);

[xrRef] = solveIncompFlow(  xrRef, [], G, S, fluid, 'bc', bc);
[xrMs]  = solveIncompFlowMS(xrMs,  [], G, CG, p, S, CS, fluid, 'bc', bc);


% Calculate norm of flux on active boundary
normFlux = norm(xrMs .faceFlux([facesBC1; facesBC2]) - ...
   xrRef.faceFlux([facesBC1; facesBC2])) / ...
   norm(xrRef.faceFlux([facesBC1; facesBC2]));

% Plot output
if doPlot,
   plot_x = @(x) plotCellData(G, x);
   plot_p = @(x) plot_x(x.cellPressure);
   plot_v = @(x) plot_x(sqrt(S.C.' * abs(x.cellFlux)));

   subplot(2,2,1)
   plot_p(xrRef); title('Pressure Fine'), set_view()
   cax = caxis; colorbar

   subplot(2,2,2)
   plot_p(xrMs); title('Pressure Coarse'), set_view()
   caxis(cax); colorbar

   subplot(2,2,3)
   plot_v(xrRef); title('Flux intensity Fine'), set_view()
   cax2 = caxis; colorbar

   subplot(2,2,4)
   plot_v(xrMs); title('Flux intensity Coarse'), set_view()
   caxis(cax2); colorbar
end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function set_view()
   view(3), axis tight, axis equal
end
