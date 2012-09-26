% SimpleModelHM_small
% a 5 x 5 reservoir for history matching
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 5; ny = 5; nz = 1;

rock.poro = repmat(0.3, [G.cells.num, 1]);

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
fluid      = adjointFluidFields(fluid);
resSolInit = initResSol(G, 0.0);

S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

totRate = 1/day;

% Introduce wells
radius = .1;

W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'rate' ,'Val',  totRate   , 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 1         , 'Type', 'rate', 'Val',   -totRate/4  , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, nx        , 'Type', 'rate', 'Val',   -totRate/4  , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'rate', 'Val',   -totRate/4  , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny     , 'Type', 'rate', 'Val',   -totRate/4  , 'Radius', radius, 'Name', 'p4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
 