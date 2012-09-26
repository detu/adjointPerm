function [W, fluid, S, resSolInit, verbose, verboseLevel] = modelCase4(G, rock)
% model Case4
% 31x41 grid blocks
%
% whether or not to show output
verbose      = false;
verboseLevel = 0;

global numControlInterval;

nx = 31; ny = 41; nz = 1;

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
fluid      = adjointFluidFields(fluid);
S          = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
W          = assembleWellSystem(G, W, 'Type', 'comp_hybrid');
resSolInit = initResSol(G, 30*barsa, 0.0);
% resSolInit = initResSol(G, 0.0);

% totVol   = sum(G.cells.volumes.*rock.poro);
% schedule = initSchedule(W, 'NumSteps', numControlInterval, 'TotalTime', 360*day, 'Verbose', verbose);
% 
% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                                   'Verbose', verbose, ...
%                                   'NumControlSteps', 1);
return;
                    
 
