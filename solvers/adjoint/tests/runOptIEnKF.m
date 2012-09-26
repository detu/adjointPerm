function permOpt = runOptIEnKF(permeability, g_mEns, g_m, asimStep, resData, cov)
% optimization of permeability in EnKF framework
% TO DO: -edit objective function !

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, param] = deal(resData{:});

% perform derivative checks
usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov};
Uinit   = permeability;
Uinit   = Uinit(:);
deriv_checkKF( Uinit, 1, usr_par);

% U   = permeability;
% %lb  = zeros(size(permeability),1);
% %ub  = ones(size(permeability),1);
% lb  = [];
% ub  = [];
% 
% A   = [];
% b   = [];
% 
% Aeq = [];
% beq = [];
% 
% % prior permeability value
% %param.m = permeability;
% priorPerm = param.m;
% 
% auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov, priorPerm};
% %options = optimset('MaxIter', 0);
% HessianFunction = @(varargin) HessTimesVecKF(varargin{:}, auxdata);
% options = optimset('HessMult', HessianFunction, 'MaxIter', 0);
% 
% tic;
% [U,fval,exitflag,output,lambda] = ktrlink(@(U)mismatchFw(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options,'knitro.opt');
% toc;
% 
% % tic;
% % [U,fval,exitflag,output,lambda] = fmincon(@(U)mismatchFw(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options);
% % toc;
% 
% permOpt = U;

return;