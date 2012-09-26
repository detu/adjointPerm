% Iterative EnKF implementation of 21x21 reservoir
% 12 wells: 6 injectors and 6 producers
% Measurement data are water-cut at each well
% TO DO: -variable for storing prior permeability
%        -update its value in every assimilation step

% cd /home/petra1a/suwartad/adjoint/adjoint-perm;
% startup;
% cd solvers/adjoint/tests/;

clear, clc;
warning off all;
opengl neverselect;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
addpath ../adjoint
addpath ../optimization
   
% % whether or not to show output
% verbose      = false;
% verboseLevel = 0;

% load 101 realizations, unit is in milli darcy
load permRealizations.mat;

% reservoir setting
nx = 21; ny = 21; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% rock properties
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

% define model function and objective function
modelFunction     = str2func('modelPermCase5');
objectiveFunction = str2func('waterCutEnKF');

[W, fluid, S, resSolInit, verbose, verboseLevel] = modelFunction(G, rock);

totVol        = sum(G.cells.volumes.*rock.poro);
schedule      = initSchedule(W, 'NumSteps', 20, 'TotalTime', 100*totVol, 'Verbose', verbose);
injMaxMin     = repmat( [1e-15 1e3]*meter^3/day, 1, 1);
prodMaxMin    = repmat( [-1e3 -1e-15]*meter^3/day, 4, 1);
controls      = initControls(schedule, 'ControllableWells', (1:5), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);

% parameters to be estimated
param.K = Kreal';
m       = ones(nx*ny*nz,1); %initial permeability assumes to be uniform
param.m = m;

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
permOpt   = zeros(numGrid,numberOfEnsemble);
priorPerm = zeros(numGrid,numberOfMeasment);

for k=1:numAsim
    
    % forecast step
    for j=1:numberOfEnsemble
        % run forward simulation for each realization
        % -------------------------------------------
        % rock properties
        m     = 1 ./ permeability(:,j);
        Kreal = ones(nx,ny);
        m     = reshape(m, [1,G.cells.num]);
        Kreal = reshape(Kreal, [1,G.cells.num]);
        K     = Kreal ./ m;
        
        save Kmodel.mat K;

        perm      = reshape(K, [1,G.cells.num]);
        perm      = reshape(perm', [1,G.cells.num]);
        rock.perm = perm'*100*milli*darcy;
        S         = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'InnerProduct', 'ip_tpf');
        
        % change schedule becomes 1 step
        schedule  = initSchedule(W, 'NumSteps', 1, 'TotalTime', 5*totVol, 'Verbose', verbose);
        
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
    
    % assimilation/analysis step USE OPTIMIZATION HERE !!!
    % ----------------------------------------------------
    
    % compute covariance Cd and Cm 
    % compute average of parameter and output
    meanPerm = mean(permeability(:,1:100),2);  
    meanFw   = mean(g_mEns,2);

    % compute matrix Y, the different of m and observation dobs
    mDif = permeability(:,1:100) - meanPerm*ones(1,numberOfEnsemble);
    gDif = g_mEns - meanFw*ones(1,numberOfEnsemble);
    
    % run optimization for each ensemble
    resData = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, param};
    for j = 1:numberOfEnsemble
        cov.Cd       = diag(gDif(:,j));
        cov.Cm       = diag(mDif(:,j));
        %cov.Cd       = eye(size(gDif(:,j),1)); %just for testing purposes
        %cov.Cm       = eye(size(mDif(:,j),1));
        permOpt(:,j) = runOptIEnKF(permeability(:,j), g_mEns(:,j), g_m, k, resData, cov);
    end
    
    % perlu update param.m from previous assimilation step !
    priorPerm = permOpt;

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

%save permCase4.mat permStep

%% DISPLAY PERMEABILITY, THE TRUTH AND MEAN OF REALIZATION FOR EACH
%% ASSIMILIATION STEP

% show truth permeability
clf;
figure(1);
load permRealizations.mat;
% rock properties
% m     = reshape(m, [1,G.cells.num]);
% Kreal = perm(:,76);    %truth permeability no. 76 from realization
% Kreal = reshape(Kreal, [1,G.cells.num]);
% K     = Kreal ./ m;
% save Kreal.mat Kreal;
% save Kmodel.mat K;
m         = 1 ./ perm(:,76);
Kreal     = ones(nx,ny);
m         = reshape(m, [1,G.cells.num]);
Kreal     = reshape(Kreal, [1,G.cells.num]);
K         = Kreal ./ m;

perm      = reshape(K, [1,G.cells.num]);
perm      = reshape(perm', [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
plotCellData(G, log10(rock.perm)); axis tight;

% show mean realizations
% figure(2);
load permCase4.mat;
for asimStep=1:9
    figure(asimStep+1);
%     Kreal     = permStep(:,asimStep); 
%     Kreal     = reshape(Kreal, [1,G.cells.num]);
%     K         = Kreal ./ m;
    
    m         = 1 ./ permStep(:,asimStep);
    Kreal     = ones(nx,ny);
    m         = reshape(m, [1,G.cells.num]);
    Kreal     = reshape(Kreal, [1,G.cells.num]);
    K         = Kreal ./ m;
    
    save Kreal.mat Kreal;
    save Kmodel.mat K;

    perm      = reshape(K, [1,G.cells.num]);
    perm      = reshape(perm', [1,G.cells.num]);
    rock.perm = perm'*100*milli*darcy;
    plotCellData(G, log10(rock.perm)); axis tight;
    % perlu subplot disini 
end
