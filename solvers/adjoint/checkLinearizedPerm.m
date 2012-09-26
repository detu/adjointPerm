function [linSimResNum, linSimRes] = checkLinearizedPerm()
% check linearized code by comparing to a perturbed control-vector
% u + epsilon*v  giving a perturbed state vector x + epsilon*w
% For small epsilon, w and v should then satify the linearized 
% equations F_x*w = -F_u*v

clear;
clc;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
addpath ../optimization
addpath ../adjoint

epsilon = 1e-8;

% whether or not to show output
verbose      = false;
verboseLevel = 0;

nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);
m  = rand(nx,ny);


Kreal = [ 1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1 ];
      
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;
perm  = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
modelHMSmall3  %update permeability

schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 100*day, 'Verbose', verbose);
controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% parameters to be estimated
param.K = reshape(Kreal,G.cells.num,1);
param.m = repmat(m',numel(schedule),1);

numW = numel(W);
u    = m;
u    = u(:);
dimU = length(u);

% v    = rand(dimU, 1) - .5;
v    = rand(dimU, 1) ;

%-------------- Initial and perurbed sim ----------------------------------
simRes  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

m         = u + epsilon*v;
m         = reshape(m, [1,G.cells.num]);
Kreal     = reshape(Kreal, [1,G.cells.num]);
K         = Kreal ./ m;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
modelHMSmall3  
simRes_p  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

linSimResNum = simRes;
for k = 2 : controls.numControlSteps+1
    linSimResNum(k).resSol.cellPressure = ...
        (simRes_p(k).resSol.cellPressure-simRes(k).resSol.cellPressure)/epsilon;
    linSimResNum(k).resSol.facePressure = ...
        (simRes_p(k).resSol.facePressure-simRes(k).resSol.facePressure)/epsilon;
    linSimResNum(k).resSol.cellFlux = ...
        (simRes_p(k).resSol.cellFlux-simRes(k).resSol.cellFlux)/epsilon;
    linSimResNum(k).resSol.faceFlux = ...
        (simRes_p(k).resSol.faceFlux-simRes(k).resSol.faceFlux)/epsilon;
    linSimResNum(k).resSol.s = ...
        (simRes_p(k).resSol.s-simRes(k).resSol.s)/epsilon;
    for k1 = 1:numW
        linSimResNum(k).wellSol(k1).flux = ...
            (simRes_p(k).wellSol(k1).flux-simRes(k).wellSol(k1).flux)/epsilon;
        linSimResNum(k).wellSol(k1).pressure = ...
            (simRes_p(k).wellSol(k1).pressure-simRes(k).wellSol(k1).pressure)/epsilon;
    end
end
    

%-------------- Linearized model -----------------------------------------
m         = u;
m         = reshape(m, [1,G.cells.num]);
Kreal     = reshape(Kreal, [1,G.cells.num]);
K         = Kreal ./ m;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
modelHMSmall3  

deltaU    = v;
linSimRes = runLinearizedPerm1(simRes, G, S, W, rock, fluid, schedule, deltaU,  'VerboseLevel', 0);

%-------------- Check output if nargout == 0 ---------------------------

if nargout == 0
    disp(['Relative differences (epsilon = ' num2str(epsilon) ' ):'])
    for k = 2 : controls.numControlSteps+1
        ep = norm(linSimResNum(k).resSol.cellPressure - linSimRes(k).resSol.cellPressure)./...
            norm(linSimResNum(k).resSol.cellPressure);
        ef = norm(linSimResNum(k).resSol.cellFlux - linSimRes(k).resSol.cellFlux)./...
            norm(linSimResNum(k).resSol.cellFlux);
        es = norm(linSimResNum(k).resSol.s - linSimRes(k).resSol.s)./...
            norm(linSimResNum(k).resSol.s);
        disp(['Step ' num2str(k-1) ' : ' num2str(ep) '  ' num2str(ef) '  ' num2str(es)])
    end
end