% function [ ]  = deriv_check( x, hess, usr_par)
% Check derivatives using finite differences

function [ ]  = deriv_checkKF( x, hess, usr_par)

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov] = deal(usr_par{:});

% FORWARD SOLVE
priorPerm = param.m;
param.m   = x;
simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);

% % UPDATE g_mEns HERE !!!
% % collected water cut as measurement
% resSol    = simRes(asimStep+1).resSol;
% wellSol   = simRes(asimStep+1).wellSol;
% [wellRates, rateSigns] = getRates(W, wellSol);
% wellCells = vertcat( W.cells );
% wellSats  = resSol.s( wellCells );
% 
% f_w_all   = fluid.fw(resSol);
% f_w       = f_w_all(wellCells);
% f_o       = 1 - f_w;
% injInx    = (rateSigns > 0);
% prodInx   = (rateSigns < 0);
% g_mEns    = f_w( prodInx );

% ADJOINT SOLVE
%adjRes    = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction);

adData    = {g_mEns, g_m, asimStep, cov, priorPerm};
adjRes    = runAdjointPermKF(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction, adData);

% COMPUTE GRADIENT
%grad = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
grad = computeGradientDKF(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, adData);
du   = cell2mat( {grad{:}}');
% du   = du';
    
fprintf(1,' Derivative checks \n')

fprintf(1,' Gradient check using finite differences (FDs)\n')
fprintf(1,' FD step size    <grad,v>     FD approx.   absolute error \n')

load Kreal;
%obj    = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule);
obj    = objectiveFunction(G, W, fluid, g_mEns, g_m, asimStep, cov, priorPerm, param);
f      = obj.val;
g      = du;
dir    = rand(size(x));
dg     = xprod(g, dir, usr_par);
delta  = 1.e-3;

for d = 1:9
    delta     = delta/10;
    uCur      = x + delta*dir;
    m         = uCur;
    param.m   = m;
    simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);
    %obj       = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule);
    
    % collected water cut as measurement
    resSol    = simRes(asimStep+1).resSol;
    wellSol   = simRes(asimStep+1).wellSol;
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells );
    
    f_w_all   = fluid.fw(resSol);
    f_w       = f_w_all(wellCells);
    f_o       = 1 - f_w;
    injInx    = (rateSigns > 0);
    prodInx   = (rateSigns < 0);
    g_mEns    = f_w( prodInx );
    obj       = objectiveFunction(G, W, fluid, g_mEns, g_m, asimStep, cov, priorPerm, param);
    f1        = obj.val;

    fprintf(1,' %12.6e  %12.6e  %12.6e  %12.6e  \n', ...
        delta, dg, (f1-f)/delta, abs(dg - (f1-f)/delta) )
    
end


if( hess )
    % perform Hessian check
    fprintf(1,'\n Hessian check using finite differences (FDs)\n')
    fprintf(1,' FD step size   absolute error \n')
    g       = du;
    dir     = rand(size(x));
    delta   = 1.e-4;
    %usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction};
    usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov, priorPerm};
    h       = HessvecKF( dir, x, usr_par );
    for d = 1:9
        delta     = delta/10;
        uCur      = x + delta*dir;
        
        m         = uCur;
        param.m   = m;
        simRes    = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);
        %adjRes    = runAdjointPerm(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction);
        
        % collected water cut as measurement
        resSol    = simRes(asimStep+1).resSol;
        wellSol   = simRes(asimStep+1).wellSol;
        [wellRates, rateSigns] = getRates(W, wellSol);
        wellCells = vertcat( W.cells );
        wellSats  = resSol.s( wellCells );

        f_w_all   = fluid.fw(resSol);
        f_w       = f_w_all(wellCells);
        f_o       = 1 - f_w;
        injInx    = (rateSigns > 0);
        prodInx   = (rateSigns < 0);
        g_mEns    = f_w( prodInx );
        
        adData    = {g_mEns, g_m, asimStep, cov, priorPerm};
        adjRes    = runAdjointPermKF(simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, modelFunction, adData);
        %grad      = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param);
        grad      = computeGradientDKF(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction, param, adData);
        g1        = cell2mat( {grad{:}}');
        
        err     = xprod(h - (g1-g)/delta, h - (g1-g)/delta, usr_par);
        aaa     = (g1-g)/delta;
        fprintf(1,' %12.6e  %12.6e   \n', delta, err )
    end

    % check selfadjointness of Hessian
    fprintf(1,'\n Check if Hessian is selfadjoint \n')

    uCur      = x + delta*dir;
    m         = uCur;
    param.m   = m;
   
    %usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction};
    usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, modelFunction, g_mEns, g_m, asimStep, param, cov, priorPerm};
    v1      = rand(size(x));
    v2      = rand(size(x));
    h       = HessvecKF( v1, x, usr_par );
    prod1   = xprod(h, v2, usr_par);
    h       = HessvecKF( v2, x, usr_par );
    prod2   = xprod(h, v1, usr_par);
    fprintf(1,'<H*v1,v2> = %12.6e, <H*v2,v1> =  %12.6e   \n', prod1, prod2 )
    fprintf(1,'| <H*v1,v2> - <H*v2,v1> | =  %12.6e   \n', abs(prod1 - prod2) )

end
fprintf(1,' End of derivative checks \n\n')
