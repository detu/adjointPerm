function simRes_refSmall = numCase4True(numControlInterval)
% True reservoir for numCase2
% a 31 x 41 reservoir for history matching
% whether or not to show output
verbose      = false;
verboseLevel = 0;


nx = 31; ny = 41; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

load permRealization;

m     = reshape(m, [1,G.cells.num]);
Kreal = perm6;
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
save Kmodel.mat K;

perm      = reshape(K, [1,G.cells.num]);
perm      = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);


% Introduce wells
radius  = .1;
totRate = 1/day;
W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'rate' ,'Val',  totRate   , 'Radius', radius, 'Name', 'i1');
W = addWell(G, rock, W, 1         , 'Type', 'bhp', 'Val',    20*barsa  , 'Radius', radius, 'Name', 'p1');
W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',    20*barsa  , 'Radius', radius, 'Name', 'p2');
W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',    20*barsa  , 'Radius', radius, 'Name', 'p3');
W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',    20*barsa  , 'Radius', radius, 'Name', 'p4');

%totRate = 4/day;
%W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'rate' ,'Val',  totRate   , 'Radius', radius, 'Name', 'i1');
%W = addWell(G, rock, W, 1         , 'Type', 'bhp', 'Val',    -0.25*totRate  , 'Radius', radius, 'Name', 'p1');
%W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',    -0.25*totRate  , 'Radius', radius, 'Name', 'p2');
%W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',    -0.25*totRate  , 'Radius', radius, 'Name', 'p3');
%W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',    -0.25*totRate  , 'Radius', radius, 'Name', 'p4');

% W = addWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'bhp' ,'Val',  35*barsa   , 'Radius', radius, 'Name', 'i1');
% W = addWell(G, rock, W, 1         , 'Type', 'bhp', 'Val',    28*barsa  , 'Radius', radius, 'Name', 'p1');
% W = addWell(G, rock, W, nx        , 'Type', 'bhp', 'Val',    28*barsa  , 'Radius', radius, 'Name', 'p2');
% W = addWell(G, rock, W, nx*ny-nx+1, 'Type', 'bhp', 'Val',    28*barsa  , 'Radius', radius, 'Name', 'p3');
% W = addWell(G, rock, W, nx*ny     , 'Type', 'bhp', 'Val',    28*barsa  , 'Radius', radius, 'Name', 'p4');
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
Wref       = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 30*barsa, 0.0);
% resSolInit = initResSol(G, 0.0);

schedule        = initSchedule(Wref, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);
simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              
