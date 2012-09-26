function [G, S, rock, W, fluid] = fiveSpot(permcase, mobcase)
% whether or not to show output

load test/data/silvia.mat;
nx = 31; ny = 41; nz = 1;
cellDims = [nx ny nz];
physDims = [nx ny nz]*2;

G = cartGrid(cellDims, physDims);   
G = computeGeometry(G);

switch permcase
    case 1
        rock.perm  = ones(G.cells.num);
        rock.poros = repmat(0.2, [G.cells.num, 1]);
    case 2
%         load rock2
%         rock.perm = 100*perm;
%         rock.poros = poros;
        rock.poros = repmat(0.2, [G.cells.num, 1]);

end

switch mobcase
    case 1
        muw = 1; muo = 1;
        totRate = 1.4;
    case 2
        muw = 5; muo = 1;
        totRate = .5;
    case 3
%         muw = 1; muo = 5;
        muw = 0.3; muo = 3;
        totRate = 10;
end

% fluid  = initTwoPhaseIncompressibleFluid('muw', muw, 'muo', muo);
fluid  = initTwoPhaseIncompressibleFluid('muw', muw, 'muo', muo, 'swc', .2, 'sor', .2);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');

% Introduce wells
radius = .1;
% W = putWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2)     , 'Type', 'bhp' ,'Val',    1 , 'Radius', radius, 'Name', 'p1');
% W = putWell(G, rock, W, 1             , 'Type', 'bhp', 'Val',    .2 , 'Radius', radius, 'Name', 'i1');
% W = putWell(G, rock, W, nx            , 'Type', 'bhp', 'Val',    .3 , 'Radius', radius, 'Name', 'i2');
% W = putWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val',    .3 , 'Radius', radius, 'Name', 'i3');
% W = putWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val',    .4 , 'Radius', radius, 'Name', 'i4');

W = putWell(G, rock,[], nx*floor(ny/2) + ceil(nx/2) +5 , 'Type', 'rate' ,'Val',  totRate , 'Radius', radius, 'Name', 'p1');
W = putWell(G, rock, W, 1             , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'i1');
W = putWell(G, rock, W, nx            , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'i2');
W = putWell(G, rock, W, nx*ny-nx+1    , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'i3');
W = putWell(G, rock, W, nx*ny         , 'Type', 'rate', 'Val',      -totRate/4 , 'Radius', radius, 'Name', 'i4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');                          
return;                              
                              
 