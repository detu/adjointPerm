function simRes_refSmall = numCase3True(numControlInterval)
% True reservoir for numCase2
% a 15 x 15 reservoir for history matching
% whether or not to show output
verbose      = false;
verboseLevel = 0;


nx = 15; ny = 15; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

load perm15;

m     = reshape(m, [1,G.cells.num]);
Kreal = rot90(perm);
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
% injectors
I   =  8 ;
J   =  8;
R   = 100/day;
nIW = 1:numel(I); 
W   = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), nz:nz, 'Type', 'rate', ...
                    'Val', R(i), 'Radius', 0.1, 'Comp_i', [1,0,0], ...
                    'name', ['I$_{', int2str(i), '}$']);
end

% producers
I   = [ 3, 3 , 12, 12 ];
J   = [ 4, 13,  4 , 13 ];
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1, 'Type', 'rate', ...
                    'Val', -25/day , 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$']);
end

fluid      = initCoreyFluid('mu', [1 5], 'sr', [0 0] );
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
Wref       = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);
% welSolInit = initWellSol(W, 30*barsa());


schedule        = initSchedule(Wref, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);

simRes_refSmall = runSchedule(resSolInit, G, S, Wref, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

save simResSmallRef.mat simRes_refSmall;
                          
return;                              