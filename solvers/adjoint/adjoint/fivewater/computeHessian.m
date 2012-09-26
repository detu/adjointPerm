function hess = computeHessian(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction)
% compute Hessian by 2.order adjoint, each column obtained by H * e_k  


%-------------- Simulation and adjoint-------------------------------------
simRes  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
adjRes  = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', 0);


%-------------- Apply hessian to each unit vector -------------------------
u    = [controls.well.values]';
u    = u(:);
dimU = length(u);
hess = zeros(dimU, dimU);

h = waitbar(0,'Computing Hessian times unit vectors ...');
for k = 1:dimU
    e_k       = zeros(dimU, 1); e_k(k) = 1;
    hess(:, k) = hessVec(G, S, W, rock, fluid, simRes, adjRes, ...
        schedule, controls, objectiveFunction, 'Vector', e_k);
    waitbar(k/dimU,h)
end
close(h)

