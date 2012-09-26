function [hessNum, hess] = checkHessian()
% check Hessian obtained by second order adjoint code by comparing to 
% numerical Hessian obtained by perturbing gradient (by adjoint) at 
% current u by epsilon*e_k , k = 1,...,dimU.  

initSimpleModel

hn  = computeNumericalHessian(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction);
h   = computeHessian(G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction);


disp(['Relative norm of differerence in Hessians: ' num2str(norm(hn-h)/norm(hn))]);
disp(['Relative norm of non-symmetric part: ' num2str(norm(h-h')/norm(h))]);

hessNum = hn;
hess    = h;