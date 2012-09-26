% test_tpf.m
clear all, clear classes

fluid = FluidProperties(SimpleCompressibleFluid());
grid  = CartGrid([4,4,1]).computeGeometry();

rock.perm=ones(grid.cells.num,1);

resSol  = initResSol (grid, 0);
wellSol = initWellSol(grid, 0);


boundary = Boundary(grid);
boundary.pside('RIGHT', (1:4), (1:1), 20);
boundary.pside('LEFT',  (1:4), (1:1), 10);

psolver = TwoPointFluxSolver(grid, rock, boundary);
[resSol, wellSol, A, rhs] = psolver.solve(resSol, wellSol, fluid);