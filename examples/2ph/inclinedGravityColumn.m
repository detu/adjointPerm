clear

%%
n      = [ 30,  30,  50];  % -> 45,000 cells
box_sz = [100, 100, 200];
G = computeGeometry(cartGrid(n, box_sz));

%% Set rock and fluid parameters
rock  = struct('perm', repmat(0.1*darcy, [G.cells.num, 1]), ...
               'poro', repmat(0.3      , [G.cells.num, 1]));

fluid = initSimpleFluid('mu', [0.307, 0.049], 'rho', [973, 617]);

%% Set gravity direction
rot   = @(theta) makehgtform('xrotate',  theta(1), ...
                             'yrotate', -theta(2), ...
                             'zrotate', -theta(3));
mul   = @(a,b,n) a(1:n,1:n) * reshape(b(1:n), [], 1);
angle = [pi/4, pi/6, 0];

gravity reset on
gravity(mul(rot(angle), gravity(), 3));

%% Assemble linsys components
S = computeMimeticIP(G, rock, 'verbose', true);

%% Solve flow problem
% Put region of CO2 at bottom of reservoir.
xr = initResSol(G, 1*barsa, 1);
d  = gravity() ./ norm(gravity);
c  = G.cells.centroids * d.' > 110;
xr.s(c) = 0;
xr = solveIncompFlow(xr, [], G, S, fluid);

%% Plot flow output
cla reset
ax = gca;
h = plotCellData(G, convertTo(xr.cellPressure, barsa), ...
                 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.625);
trans = hgtransform;
set(trans, 'Parent', ax, 'Matrix', ...
    makehgtform('translate', [ box_sz(1), 0,  box_sz(3)],    ...
                'xrotate'  , -angle(1), 'yrotate', angle(2), ...
                'translate', [-box_sz(1), 0, -box_sz(3)]));
set(h, 'Parent', trans);
view([45, 5]), grid on, camproj perspective
colorbar

mat = get(trans, 'Matrix');

%%
dT = [1, 2, 2, 5, 5, 10, 15, 20, 40, 50, 50, ...
      100, 100, 200, 200, 300, 400, 500] .* day();
dT = [dT, [2, 2, 2, 4, 5, 5, 10, 10, repmat(15, [1, 10])].*year()];

clf, ax = gca;
h = plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
trans = hgtransform;
set(trans, 'Parent', ax, 'Matrix', mat);
set(h, 'Parent', trans);

thresh = 0.999;
hs = plotCellData(G, xr.s, find(xr.s < thresh), 'EdgeColor', 'k', ...
                  'EdgeAlpha', 0.1, 'FaceAlpha', 0.25);
set(hs, 'Parent', trans)

s  = linspace(0, 1, 8192).';
cm = [1-s.^(13/16), 1-s.^6, s.^6];
colormap(cm)

view([45, 5]), grid on, camproj perspective
colorbar
%%
t = 0;
for k = 1 : numel(dT),
   xr = explicitTransport(xr, [], G, dT(k), rock, fluid, ...
                          'Verbose', true, 'SatWarn', 1.0e-3);

   % Check for inconsistent saturations
   assert (max(xr.s) < 1+eps && min(xr.s) > -eps);

   % Increase time and plot saturation
   t = t + dT(k);
   delete(hs)
   hs = plotCellData(G, xr.s, find(xr.s < thresh), 'EdgeColor', 'k', ...
                     'EdgeAlpha', 0.1, 'FaceAlpha', 0.25);
   set(hs, 'Parent', trans)

   drawnow

   % Compute new flow field.
   xr = solveIncompFlow(xr, [], G, S, fluid);
end
