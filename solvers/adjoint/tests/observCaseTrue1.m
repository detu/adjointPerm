function simRes_refSmall = observCaseTrue1(numControlInterval)
% True reservoir for observCase
% a 17 x 17 reservoir for history matching
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 17; ny = 17; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

load spe10.mat
perm = KU(1,1:nx,1:ny,9);
perm = reshape(perm,[G.cells.num, 1]);

m     = reshape(m, [1,G.cells.num]);
Kreal = perm / (100*milli*darcy);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
save Kmodel.mat K;

perm = reshape(K, [1,G.cells.num]);

perm = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

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
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
Wref       = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);

schedule        = initSchedule(Wref, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);

simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              