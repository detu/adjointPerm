function [W, fluid, S, resSolInit, verbose, verboseLevel] = modelObs(G, rock)
% a 25 x 1 reservoir for history matching

% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 25; 
ny = 1; 
nz = 1;

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
fluid      = adjointFluidFields(fluid);
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
W          = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 300*barsa);

return;
                    
 