function hess = computeHessianPerm(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param)
% compute Hessian by 2.order adjoint, each column obtained by H * e_k  


%-------------- Simulation and adjoint-------------------------------------
simRes  = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, 'VerboseLevel', 0);
adjRes  = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param,  'VerboseLevel', 0);


%-------------- Apply hessian to each unit vector -------------------------
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
hess     = zeros(dimU, dimU);

h = waitbar(0,'Computing Hessian times unit vectors ...');
for k = 1:dimU
    e_k        = zeros(dimU, 1); e_k(k) = 1;
    hess(:, k) = hessVecPerm(G, S, W, rock, fluid, simRes, adjRes, schedule, controls, objectiveFunction, param, 'Vector', e_k);
    waitbar(k/dimU,h)
end
close(h)

