% EnKF implementation of 45x45 reservoir
% 12 wells: 6 injectors and 6 producers
% Measurement data are water-cut at each well


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
load permeability;

% reservoir setting
nx = 45; ny = 45; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = ones(nx,ny);

% rock properties
m     = reshape(m, [1,G.cells.num]);
Kreal = permeability(:,101);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

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
W = [];
W = addWell(G, rock,W, 1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', 'I1');
for i=1:5
    W = addWell(G, rock,W, (9*i - 1)*nx +1      , 'Type', 'rate' ,'Val',  totRate  ,'Radius', radius, 'Name', ['I',num2str(i+1)]);
end

W = addWell(G, rock,W, 45     , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', 'P1');
for i=1:5
    W = addWell(G, rock,W, i*9*nx      , 'Type', 'rate' ,'Val',  -totRate  ,'Radius', radius, 'Name', ['P',num2str(i+1)]);
end

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol     = sum(G.cells.volumes.*rock.poro);
schedule   = initSchedule(W, 'TimeSteps', 120*day:120*day:1200*day, 'Verbose', verbose);

injMaxMin  = repmat( [0.001 100]*meter^3/day, 6, 1);
prodMaxMin = repmat( [-100 -0.001]*meter^3/day, 6, 1);

controls   = initControls(schedule, 'ControllableWells', (1:12), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', [], ...
                       'Verbose', verbose);
load Uopt;
controls  = updateControls(controls, Uopt);
schedule  = updateSchedule(controls, schedule);
                   
%% -----------------------------------------------------------------------------------------
%  EnKF part
g_m         = zeros(12,10);
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
    g_m(:,step-1) = f_w;
end

% number of measurement
numAsim   = 9;
g_mEns    = zeros(12,100);
numGrid   = nx*ny*nz;
numberOfEnsemble = 100;
numberOfMeasment = 12;
mG        = zeros(2025,12);
mGEns     = zeros(2025,12);
gG        = zeros(12,12);
gGEns     = zeros(12,12);
permStep  = zeros(numGrid,numAsim);
cD        = 0.005;

for k=1:numAsim
    
    % forecast step
    for j=1:numberOfEnsemble
        % run forward simulation for each realization
        % -------------------------------------------
        % rock properties
        m     = ones(nx,ny);
        m     = reshape(m, [1,G.cells.num]);
        Kreal = permeability(:,j);
        Kreal = reshape(Kreal, [1,G.cells.num]);
        K     = Kreal ./ m;

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
        g_mEns(:,j) = f_w;
        
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
        yUpd   = [permeability(:,j)*100*milli*darcy ; g_mEns(:,j)] + ... 
                  Ke*(g_m(:,k) + vn*ones(numberOfMeasment,1)  - g_mEns(:,j));

        % update permeability ensemble
        permeability(:,j) = yUpd(1:numGrid) / 100*milli*darcy;

    end

    % store average permeability for display/checking purpose
    permStep(:,k)  = mean(permeability(:,1:100),2);
    

    
end

%% DISPLAY TRUTH PERMEABILITY AND ITS ESTIMATED (ENSEMBLE AVERAGE) FOR EACH ASSIMILATION STEP