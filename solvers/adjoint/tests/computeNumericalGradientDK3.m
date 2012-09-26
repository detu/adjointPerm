function numGrad = computeNumericalGradientDK3(simRes, G, S, W, rock, fluid, schedule, objectiveFunction, param, resSolInit, modelFunction)
% compute numerical gradient

epsilon   = 1e-6;
obj       = objectiveFunction(param, G, S, W, rock, fluid, simRes);
valInit   = obj.val;
nx        = G.cartDims(1); 
ny        = G.cartDims(2); 
nz        = G.cartDims(3);
gd        = nx * ny* nz;

load Kreal;
m        = param.m;
uInit    = m(:);
dimU     = length( uInit );
values   = zeros(1,dimU);

for k = 1:dimU
    e_k       = zeros(dimU, 1); 
    e_k(k)    = 1;
    uCur      = uInit + epsilon*e_k;
    param.m   = uCur ;
%     simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, 'VerboseLevel', 0);
    simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);
    obj       = objectiveFunction(param, G, S, W, rock, fluid, simRes);
    values(k) = obj.val;
end

numGrad = reshape((values' - valInit)/epsilon, 1, gd);

  