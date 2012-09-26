function simRes_refSmall = observCaseTrue(numControlInterval)
% True reservoir for observCase
% a 25 x 1 reservoir for history matching
% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 25; ny = 1; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

load perm;
m     = reshape(m, [1,G.cells.num]);
Kreal = perm;
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
% W = addWell(G, rock,[], 1         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd1');
% W = addWell(G, rock, W, 13        , 'Type', 'bhp',  'Val',  320*barsa, 'Radius', radius, 'Name', 'Inj1');
% W = addWell(G, rock, W, nx        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd2');
% -- Well Configuration 2
% W = addWell(G, rock,[], 1         , 'Type', 'bhp',  'Val',  320*barsa, 'Radius', radius, 'Name', 'Inj1');
% W = addWell(G, rock, W, nx        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd1');
% -- Well Configuration 3
% W = addWell(G, rock,[], 1         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd1');
% W = addWell(G, rock, W, 5         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd3');
% W = addWell(G, rock, W, 13        , 'Type', 'bhp',  'Val',  320*barsa, 'Radius', radius, 'Name', 'Inj1');
% W = addWell(G, rock, W, 20        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd4');
% W = addWell(G, rock, W, nx        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd2');
% -- Well Configuration 4
W = addWell(G, rock,[], 1         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd1');
W = addWell(G, rock, W, 5         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd3');
W = addWell(G, rock, W, 8         , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd5');
W = addWell(G, rock, W, 11        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd7');
W = addWell(G, rock, W, 13        , 'Type', 'bhp',  'Val',  320*barsa, 'Radius', radius, 'Name', 'Inj1');
W = addWell(G, rock, W, 15        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd8');
W = addWell(G, rock, W, 18        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd6');
W = addWell(G, rock, W, 22        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd4');
W = addWell(G, rock, W, nx        , 'Type', 'bhp',  'Val',  300*barsa, 'Radius', radius, 'Name', 'Prd2');

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
Wref       = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);

schedule        = initSchedule(Wref, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);

simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              