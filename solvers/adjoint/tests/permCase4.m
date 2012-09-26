% EnKF implementation of 45x45 reservoir
% 12 wells: 6 injectors and 6 producers
% Measurement data are water-cut at each well

cd /home/petra1a/suwartad/adjoint/adjoint-perm;
startup;
cd solvers/adjoint/tests/;

clear, clc;
warning off all;
opengl neverselect;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
% whether or not to show output
verbose = false;
verboseLevel = 0;

% load 101 realizations, unit is in milli darcy
load permRealizations.mat;

% reservoir setting
nx = 21; ny = 21; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
%m  = ones(nx,ny);

% rock properties
%m     = reshape(m, [1,G.cells.num]);
%Kreal = perm(:,76);    %truth permeability no. 76 from realization
%Kreal = reshape(Kreal, [1,G.cells.num]);
%K     = Kreal ./ m;
m         = 1 ./ perm(:,76);
Kreal     = ones(nx,ny);
m         = reshape(m, [1,G.cells.num]);
Kreal     = reshape(Kreal, [1,G.cells.num]);
K         = Kreal ./ m;

% change permeability realizations
permeability       = perm;
permeability(:,76) = [];

save Kreal.mat Kreal;
save Kmodel.mat K;

perm      = reshape(K, [1,G.cells.num]);
perm      = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

% fluid properties
fluid         = initCoreyFluid('mu', [1, 5], 'sr', [0.1 0.1]);
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

totVol    = sum(G.cells.volumes.*rock.poro);
schedule  = initSchedule(W, 'NumSteps', 20, 'TotalTime', 100*totVol, 'Verbose', verbose);

injMaxMin  = repmat( [1e-15 1e3]*meter^3/day, 1, 1);
prodMaxMin = repmat( [-1e3 -1e-15]*meter^3/day, 4, 1);

controls   = initControls(schedule, 'ControllableWells', (1:5), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);

% load optimized controls
load UoptCase4.mat;
controls  = updateControls(controls, Uopt);
schedule  = updateSchedule(controls, schedule);

%% -----------------------------------------------------------------------------------------
%  EnKF part
g_m         = zeros(4,20);
% run the truth model
simResTruth = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
numSteps    = numel(simResTruth);

for step=2:numSteps
    % collected water cut as measurement
    resSol    = simResTruth(step).resSol;
    wellSol   = simResTruth(step).wellSol;
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells );
    
    f_w_all   = fluid.fw(resSol);
    f_w       = f_w_all(wellCells);
    f_o       = 1 - f_w;
    injInx    = (rateSigns > 0);
    prodInx   = (rateSigns < 0);
    g_m(:,step-1) = f_w( prodInx );
end

% number of measurement
numAsim   = 20;
numGrid   = nx*ny*nz;
numberOfEnsemble = 100;
numberOfMeasment = 4;
g_mEns    = zeros(numberOfMeasment,numberOfEnsemble);
mG        = zeros(numGrid,numberOfMeasment);
mGEns     = zeros(numGrid,numberOfMeasment);
gG        = zeros(numberOfMeasment,numberOfMeasment);
gGEns     = zeros(numberOfMeasment,numberOfMeasment);
permStep  = zeros(numGrid,numAsim);
cD        = 0.005;

for k=1:numAsim

    fprintf(1, 'Assimilation step :%d \n', k);
    
    % forecast step
    for j=1:numberOfEnsemble
        % run forward simulation for each realization
        % -------------------------------------------
        % rock properties
        %m     = ones(nx,ny);
        %m     = reshape(m, [1,G.cells.num]);
        %Kreal = permeability(:,j);
        %Kreal = reshape(Kreal, [1,G.cells.num]);
        %K     = Kreal ./ m;

        m         = 1 ./ permeability(:,j);
        Kreal     = ones(nx,ny);
        m         = reshape(m, [1,G.cells.num]);
        Kreal     = reshape(Kreal, [1,G.cells.num]);
        K         = Kreal ./ m;

        save Kreal.mat Kreal;
        save Kmodel.mat K;

        perm      = reshape(K, [1,G.cells.num]);
        perm      = reshape(perm', [1,G.cells.num]);
        rock.perm = perm'*100*milli*darcy;
        S         = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');
        simResEns = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
        step      = k + 1;
        
        % collected water cut as measurement
        resSol    = simResEns(step).resSol;
        wellSol   = simResEns(step).wellSol;
        [wellRates, rateSigns] = getRates(W, wellSol);
        wellCells = vertcat( W.cells );
        wellSats  = resSol.s( wellCells );
        
        f_w_all   = fluid.fw(resSol);
        f_w       = f_w_all(wellCells);
        f_o       = 1 - f_w;
        injInx    = (rateSigns > 0);
        prodInx   = (rateSigns < 0);
        g_mEns(:,j) = f_w( prodInx );
        
    end
    
    % assimilation/analysis step
    % -------------------------------

    % compute average of parameter and output
    meanPerm = mean(permeability(:,1:100),2);    % make sure later if it is correct
    meanFw   = mean(g_mEns,2);

    % compute matrix Y, the different of m and observation dobs
    mDif = permeability(:,1:100) - meanPerm*ones(1,numberOfEnsemble);
    gDif = g_mEns - meanFw*ones(1,numberOfEnsemble);

    % compute Kalman gain, defined as Ke
    for j = 1:numberOfEnsemble
        mG = mDif(:,j) * gDif(:,j)';
        gG = gDif(:,j) * gDif(:,j)';

        % summation for each ensemble
        mGEns = mG + mGEns;
        gGEns = gG + gGEns;
    end
    dYA = [mGEns;gGEns] / (numberOfEnsemble - 1);
    AtA = ( gGEns / (numberOfEnsemble - 1) ) + cD*eye(numberOfMeasment,numberOfMeasment);
    Ke  = dYA / AtA;

    % update state vector
    for j = 1:numberOfEnsemble

        % perturbed measurement noise
        vn = normrnd(0,cD);

        % ensemble model
        yUpd   = [permeability(:,j) ; g_mEns(:,j)] + ...
                  Ke*(g_m(:,k) + vn*ones(numberOfMeasment,1)  - g_mEns(:,j));

        % update permeability ensemble
        permeability(:,j) = yUpd(1:numGrid);

    end

    % store average permeability for display/checking purpose
    permStep(:,k)  = mean(permeability(:,1:100),2);
    
end

save permCase4.mat permStep
