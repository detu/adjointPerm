function [W, fluid, S, resSolInit, verbose, verboseLevel] = modelObs1(G, rock)
% a 17 x 17 reservoir for history matching

% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 17; 
ny = 17; 
nz = 1;

% Introduce wells
radius = .1;
% -- Well Configuration 1
% W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)    , 'Type', 'bhp' ,'Val',  310*barsa ,   'Radius', radius, 'Name', 'i1');
% W = addWell(G, rock, W, 1             , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p1');
% W = addWell(G, rock, W, nx            , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p2');
% W = addWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p3');
% W = addWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p4');
% -- Well Configuration 2
W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)    , 'Type', 'bhp' ,'Val',  300*barsa ,   'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, 1             , 'Type', 'bhp', 'Val',      310*barsa , 'Radius', radius, 'Name', '11');
W = addWell(G, rock, W, nx            , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val',      300*barsa , 'Radius', radius, 'Name', 'p4');

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
fluid      = adjointFluidFields(fluid);
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
W          = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);

return;
                    
 