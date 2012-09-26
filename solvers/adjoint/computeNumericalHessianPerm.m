function hess = computeNumericalHessianPerm(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param)
% compute numerical Hessian by comparing gradient (by adjoint) at current u
% and u+epsilon*e_k , k = 1,...,dimU 

epsilon = 1e-8;

%-------------- Initial gradient ------------------------------------------
simRes  = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, 'VerboseLevel', 0);
adjRes  = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, param, 'VerboseLevel', 0);
grad    = computeGradientPerm(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
gradVec = cell2mat( {grad{:}}');

%-------------- Gradients for perturbed controls --------------------------

% u    = [controls.well.values]';
% u    = u(:);
% dimU = length(u);

nx = G.cartDims(1); 
ny = G.cartDims(2); 
nz = G.cartDims(3);
gd = nx * ny* nz;
numSteps  = numel( schedule);

load Kreal;
m        = reshape(param.m, G.cells.num, numel(schedule));
mOrig    = m;
uInit    = m(:);
dimU     = length( uInit );
uInit    = reshape(uInit,gd,numSteps);

pertGradVecs = zeros(dimU, dimU);

h = waitbar(0,'Computing gradients for perturbed controls ...');
for k = 1:dimU
    e_k       = zeros(dimU, 1); e_k(k) = 1;
    % update permeability -> update well & assemblemimetic !
    j         = ceil(k/gd);
    Un        = uInit(:,j);
    e_k       = reshape(e_k,gd,numSteps);
    e_k       = e_k(:,j);
    uCur      = Un + epsilon*e_k;
    m         = mOrig;
    m(:,j)    = uCur;
    param.m   = m;
    Kreal     = reshape(Kreal, [G.cells.num,1]);
    rock.perm = (Kreal ./ m(:,j))*100*milli*darcy;  
    modelHMSmall3
    simRes    = runSchedulePerm(simRes(1).resSol, G, S, W, rock, fluid, schedule, param, 'VerboseLevel', 0);
    adjRes    = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, 'VerboseLevel', 0);
    grad      = computeGradientPerm(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
    pertGradVecs(:, k) = cell2mat( {grad{:}}');
    waitbar(k/dimU,h)
end
close(h)

hess = ( pertGradVecs - gradVec*ones(1,dimU) )/epsilon;


function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end    