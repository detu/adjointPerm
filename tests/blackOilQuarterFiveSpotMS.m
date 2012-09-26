% blackOilQuarterFiveSpot
%
% implisitt løser?
% blandbar strøm
% rapportér massebalanse.
% OIP,WIP,GIP,etc
%

clear

%pvtfile    = 'pvt.txt';
pvtfile    = 'compressible.txt';
%pvtfile    = 'incompressible.txt';
verbose    = true;

PVTTAB     = readpvt(fullfile(ROOTDIR, 'params', 'fluid', 'Data', pvtfile));


cellDims   = [ 40,  40,  1];
physDims   = [400, 400, 10];
G          = cartGrid(cellDims, physDims);
G          = computeGeometry(G);

rock.perm  = repmat(100.0*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(  1.0,             [G.cells.num, 1]);
%{
rock.s_w   = repmat(  0.0, [G.cells.num, 1]);
rock.s_g   = repmat(  0.0, [G.cells.num, 1]);
rock.s_o   = repmat(  1.0, [G.cells.num, 1]);
%}

S          = computeMimeticIP(G, rock);

p  = partitionUI(G, [20, 20, 1]);
p  = processPartition  (G, p, 'Verbose', verbose);
CG = generateCoarseGrid(G, p, 'Verbose', verbose);
CS = generateCoarseSystem(G, rock, S, CG, ...
                          'Verbose', verbose);

%% Wells
W   = addWell(G, rock, [], 1, 'Type', 'rate', ...
              'Val', 1000*meter^3/day, 'Radius', .1, 'Dir', 'z');
W   = addWell(G, rock, W, G.cells.num, 'Type','bhp', ...
              'Val', 420*barsa, 'Radius', .1);
W   = assembleWellSystem(G, W);
W = generateCoarseWellSystem(G, S, CG, CS, ones([G.cells.num, 1]), rock, W);

z.a = zeros([G.cells.num, 1]);
z.v = zeros([G.cells.num, 1]);
z.l = ones ([G.cells.num, 1]);

%DT = 10;
DT  = 2;

initP             = 420; %bar
resSol  = initResSol (G, initP);
wellSol = initWellSol(W, initP);

pressure_prev     = resSol.cellPressure;

fluid             = pvt(PVTTAB, rock,        ...
                        resSol.cellPressure, ...
                        S.constants.gravity);

z.a = zeros([G.cells.num, 1]);
z.v = zeros([G.cells.num, 1]);
z.l = ones ([G.cells.num, 1]) ./ fluid.Bl;
zum = fluid.Bl.*z.l+fluid.Ba.*z.a+fluid.Bv.*z.v;
volume_discrepancy = (zum - 1) .* poreVolume(G, rock) ./ DT;
S.RHS.volume_discrepancy = volume_discrepancy;
CS.RHS.volume_discrepancy = accumarray(p, volume_discrepancy);

rock.s_w = z.a;
rock.s_o = z.l;
rock.s_g = z.v;

% Splittesteg

T  = 0;
while true
  if T > 100, DT = 100; end

  %
  % Solve pressure
  %
  resV = 1;
  while resV > 0.5e-5,
    cf0 = resSol.cellFlux;

    [resSol, wellSol] = solveBlackOilWellSystemMS(resSol, wellSol,   ...
                                                  G, CG, rock, S, CS, W, ...
                                                  PVTTAB, pressure_prev, ...
                                                  DT, 'Verbose', false);

    resV = norm(resSol.cellFlux - cf0, inf);
    dispif(verbose, 'Residual in velocity: %5.5e\n', resV);
  end
  fluid = pvt(PVTTAB, rock, resSol.cellPressure, S.constants.gravity);


  % Solve saturation equation
  [gm, sources, porevol] = inflow(G, rock, S, W, resSol, wellSol);
  gp   = spdiags(-sum(gm,1)', 0, G.cells.num, G.cells.num);
  Fout = full(gp(1 : G.cells.num+1 : end)) .';
  F    = gp + gm;

  % CFL condition...
  t    = 0;
  dt   = min(0.5 .* porevol ./ abs(Fout));


  nc = 0;
  while t < DT,
    % Explicit Euler

    % Compute saturations
    zum = fluid.Bl.*z.l + fluid.Ba.*z.a + fluid.Bv.*z.v;
    s.a = fluid.Ba .* z.a ./ zum;
    s.l = fluid.Bl .* z.l ./ zum;
    s.v = fluid.Bv .* z.v ./ zum;

    kr    = getRelPerm([s.a, s.l, s.v], PVTTAB);
    fa    = kr.a ./ (kr.a + kr.l + kr.v) ./ fluid.Ba;
    fl    = kr.l ./ (kr.a + kr.l + kr.v) ./ fluid.Bl;

    % Phase flux matrices
    FA    = F*fa;
    FL    = F*fl;


    z.a = z.a - dt./porevol.*(FA - max(sources,0)*1./fluid.Ba - min(sources,0).*fa);
    z.l = z.l - dt./porevol.*(FL                              - min(sources,0).*fl);


    % Plot pressure and density
    subplot(3,1,1)
    plot(1:G.cartDims(1), resSol.cellPressure(1 : G.cartDims(1)+1 : end));
    subplot(3,1,2);
    title('z');
    h = plotyy(1:G.cartDims(1), z.a(1:G.cartDims(1)+1:end), ...
               1:G.cartDims(1), z.l(1:G.cartDims(1)+1:end)); axis(h, 'tight')

    satsum = fluid.Bl.*z.l + fluid.Ba.*z.a + fluid.Bv.*z.v;
    s.a = fluid.Ba .* z.a ./ satsum;
    s.l = fluid.Bl .* z.l ./ satsum;
    s.v = fluid.Bv .* z.v ./ satsum;

    if false
      subplot(3,1,3);
      title('s');
      h = plotyy(1:G.cartDims(1), s.a(1:G.cartDims(1)+1:end), ...
                 1:G.cartDims(1), s.l(1:G.cartDims(1)+1:end)); axis(h, 'tight')
    else
      subplot(3,1,3);
      hold off;
      surf(reshape(s.a, G.cartDims(1:2))); view(45, 55);
    end
    drawnow

    t  = t + dt;
    dt = min(dt, DT-t);
    dispif(verbose, repmat('\b', [1, nc]));
    nc = dispif(verbose, 't = %5.5f', T + t);
  end
  dispif(verbose, '\n');
  T = T + DT;

  rock.s_w = z.a;
  rock.s_g = z.v;
  rock.s_o = z.l;

  pressure_prev = resSol.cellPressure;

  zum = fluid.Bl.*z.l + fluid.Ba.*z.a + fluid.Bv.*z.v;
  volume_discrepancy = (zum - 1) .* poreVolume(G, rock) ./ DT;
  S.RHS.volume_discrepancy  = volume_discrepancy;
  CS.RHS.volume_discrepancy = accumarray(p, volume_discrepancy);
end
