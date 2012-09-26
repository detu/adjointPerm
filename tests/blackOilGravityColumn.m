clear

%%
n      = [  30,  30,  50];
box_sz = [ 100, 100, 200];
G = computeGeometry(cartGrid(n, box_sz));

%% Set rock and fluid parameters
rock  = struct('perm', repmat(0.1*darcy, [G.cells.num, 1]), ...
               'poro', repmat(0.3      , [G.cells.num, 1]));

% Initial reservoir pressure
initP =  335 * barsa();
DT    =   10 * day()  ;

fluid = initCompressibleFluid(1000, initP, 100/initP, 1*centi*poise());
s0    = 1;
z0    = 1;

%% Set gravity direction
rot   = @(theta) makehgtform('xrotate',  theta(1), ...
                             'yrotate', -theta(2), ...
                             'zrotate', -theta(3));
mul   = @(a,b,n) a(1:n,1:n) * reshape(b(1:n), [], 1);
angle = [pi/4, pi/6, 0] .* 1;

gravity reset on
gravity(mul(rot(angle), gravity(), 3));

%% Assemble linsys components
S = computeMimeticIP(G, rock, 'verbose', true);

resSol = initResSol(G, initP, s0, z0);

porvol = poreVolume(G, rock);
vd     = @(dt, u) ((sum(u, 2) - 1) .* porvol ./ dt);

[u, u, u, u] = fluid.pvt(initP, z0);

% Measure initial volume discrepancy term (== 0 by definition).
S.RHS.volume_discrepancy = vd(DT, u);

solver_press = @(xr, p0, dt) ...
   solveBlackOilWellSystem(xr, [], G, rock, S, fluid, p0, dt);

% Solve pressure equation.
%
solve_p = @(xr, varargin) solver_press(xr, resSol.cellPressure, DT);
resSol  = succSubst(resSol, [], solve_p, 'Tol', 5.0e-4/day(), ...
                    'MaxIt', 15, 'Verbose', true);
[u, u, u, u] = fluid.pvt(resSol.cellPressure, resSol.z);
resSol.s = bsxfun(@rdivide, u, sum(u, 2));

%% Plot flow output
cla reset
ax = gca;
h = plotCellData(G, convertTo(resSol.cellPressure, barsa), ...
                 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.625);
trans = hgtransform;
set(trans, 'Parent', ax, 'Matrix', ...
    makehgtform('translate', [ box_sz(1), 0,  box_sz(3)],    ...
                'xrotate'  , -angle(1), 'yrotate', angle(2), ...
                'translate', [-box_sz(1), 0, -box_sz(3)]));
set(h, 'Parent', trans);
view([45, 5]), grid on, camproj perspective
colorbar


%% Example last changed:
% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $
