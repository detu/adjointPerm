% mimetic_mixed_tpf
%
% Same as simpleCartGridExample but uses a tpf-discretization solved using
% mimetic, mixed and cell-centered systems.

% Construct 20-by-20-by-5 grid of dimensions 20-by-20-by-5 (meters)
G = cartGrid([20, 20, 5]);

% Compute areas, centroids, normals, etc ...
G = computeGeometry(G, 'Verbose', true);

rock.poro = repmat(0.3, [G.cells.num, 1]);
rock.perm  = repmat(100 * (darcy() / 1000), [G.cells.num, 1]);

fluid = initSimpleFluid();

sol_h = initResSol(G, 0.0);
sol_m = initResSol(G, 0.0);
sol_t = initResSol(G, 0.0);

xw    = initWellSol([], 0.0);

% Introduce sources and pressure BC
someCells = (210 : 20*20 : 20*20*5) .';
src = addSource([], someCells, ones(size(someCells)) ./ 86400, ...
                'sat', repmat([1,0,0], [numel(someCells), 1]));

bc  = pside([], G, 'LEFT', 1:G.cartDims(2), 1:G.cartDims(3), 1e5);

% Assemble hybrid and mixed system
Sh = computeMimeticIP(G, rock,                  ...
                           'Verbose',      true,     ...
                           'InnerProduct', 'ip_tpf');
Sm = computeMimeticIP(G, rock,                  ...
                           'Verbose',      true,     ...
                           'InnerProduct', 'ip_tpf', ...
                           'Type',         'mixed');

% Solve systems
disp('MIMETIC:')
sol_h = solveIncompFlow(sol_h, xw, G, Sh, fluid, ...
                        'bc', bc, 'src', src,    ...
                        'MatrixOutput', true);

disp('MIXED:')
sol_m = solveIncompFlow(sol_m, xw, G, Sm, fluid, ...
                        'bc', bc, 'src', src,    ...
                        'Solver', 'mixed',       ...
                        'MatrixOutput', true);
disp('TWO-POINT:')
sol_t = solveIncompFlow(sol_t, xw, G, Sm, fluid, ...
                        'bc', bc, 'src', src,    ...
                        'Solver', 'tpfa',        ...
                        'MatrixOutput', true);

% Plot output
ce = @(s  ) num2str(round(condest(s.A)));
d  = @(s,f) abs(s(1).(f) - s(2).(f));
e  = @(s,f) [d(s([1,2]), f), d(s([1,3]), f)];

figure
subplot(2,3,1)
   plot(1 : G.faces.num, e([sol_h, sol_m, sol_t], 'faceFlux') .* 86400);
   legend('Mixed', 'TPF')
   title('Flux compared to hybrid')

subplot(2,3,2)
   plot(1 : G.cells.num, e([sol_h, sol_m, sol_t], 'cellPressure') ./ 1e5);
   legend('Mixed', 'TPF')
   title('Pressure compared to hybrid')

subplot(2,3,4)
   spy( sol_h.A );
   title(['System matrix, hybrid - condest: ',       ce(sol_h)]);

subplot(2,3,5)
   spy( sol_m.A );
   title(['System matrix, mixed - condest: ',        ce(sol_m)]);

subplot(2,3,6)
   spy( sol_t.A );
   title(['System matrix, mixed (TPFA) - condest: ', ce(sol_t)]);
