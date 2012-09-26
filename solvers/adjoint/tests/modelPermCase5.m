function [W, fluid, S, resSolInit, verbose, verboseLevel] = modelPermCase5(G, rock)
% model for permCase5
% 21x21 grid blocks
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 21; ny = 21; nz = 1;

% fluid properties
fluid         = initCoreyFluid('mu', [1, 5], 'sr', [0.1 0.1]);
fluid         = adjointFluidFields(fluid);
resSolInit    = initResSol(G, 0.0);
resSolInit.s  = 0.2*ones(G.cells.num, 1);
S             = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');

% Introduce wells
radius  = .1;
totRate = 1/day;
W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'rate' ,'Val',   totRate  , 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 5*nx+6        , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, 6*nx-5        , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, 16*nx+5       , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, 17*nx-5       , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'p4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

return;
                    
 