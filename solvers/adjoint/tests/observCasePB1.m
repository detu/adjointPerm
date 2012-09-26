% Pure Barrier Method
% Analyzing of Observability 
% The analysis is similar to "output constraint" 
% see Wei Kang's paper for more detail..
% July, 15th 2010
% suwartad@itk.ntnu.no

clear, clc;
warning off all;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

addpath ../adjoint
addpath optimization

% run the TRUE model
numCtrlIntv = 1;
observCaseTrue1(numCtrlIntv);

nx = 17; ny = 17; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
load minit.mat
Kreal = ones(17,17);

m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

save Kreal.mat Kreal;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;
rock.poro = repmat(0.3, [G.cells.num, 1]);

global numControlInterval;
numControlInterval = 1;
modelFunction      = str2func('modelObs1');
[W, fluid, S, resSolInit, verbose, verboseLevel] = modelFunction(G, rock);

% Choose objective function
objectiveFunction = str2func('PermMismatchPB');
constrainFunction = str2func('Jac_Fwn');

totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 360*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:2), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);
% parameters to be estimated
param.K = Kreal';   
param.m = m';

% -- OPTIMIZATION PART --
% initial permeability value
U   = param.m;
lb  = 1e-5*ones(289,1);
ub  = 1*ones(289,1);

lbc = [];
ubc = [];

A   = []; 
b   = [];

Aeq = [];
beq = [];
 
%%% AUGMENTED BARRIER FUNCTION PARAMETERS %%%%%
% STEP 0: Initialization
% Initial Augmentend-Lagrange multipier
% should be checked that c_i(x) + s_k,i > 0 is HOLD !!!
param.tau         = 0.1;
param.mu          = 1e12;    % initial penalty parameter
% nonlinear inequality constraint
% h_bound           = 1e-2;   % output constraint with epsilon 1e-1
% initial convergence/check stopping criteria
omega             = 1e-1;   % stopping criteria for inner iteration related to norm of gradient
tau               = 1e-1;   % stopping criteria for norm of constraint violation
% final convergence/check stopping criteria
omegaStar         = 1e-8;
tauStar           = 9.5e-5;   % defined as the epsilon equals 1e-2, should be 1e-4
% initial grad of inner iteration
gradP             = inf ;

% initial (additional) convergence test (based on obj. value changes)
jk         = inf;
funcCount  = 0;   % number of forward simulation
constraint = -inf;

% TAMBAHAN STOPPING CRITERIA : constraint = 1e-4 STOP !!!
% CHECK STOPPING CRITERIA YG LAIN !

while ( gradP >= omegaStar )
    
    % STEP 1
    % inner iteration
    options = optimset('Display', 'iter', 'Hessian', 'sr1', ...
                      'TolFun', omega, 'Algorithm', 'active-set', 'GradObj', 'on', 'TolX', 1e-15, 'maxit', 50);
    auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction, param, modelFunction};
    
    [U,fval,exitflag,output] = ktrlink(@(U)ambiguityValue(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options);
    gradP     = output.firstorderopt;
    j0        = fval;
    funcCount = output.funcCount + funcCount;

    % STEP 2
    % convergence test: check constraint violation
    param.m    = U;
    simRes     = runSchedulePerm(resSolInit, G, S, W, rock, fluid, schedule, param, modelFunction, 'VerboseLevel', 0);
    obj        = constrainFunction(param, G, S, W, rock, fluid, simRes, schedule);
    constraint = obj.val;

    % - update barrier parameter
    param.mu = param.tau * param.mu;
    omega    = param.tau * omega;

    % test convergence based on obj. func. changes
    if (abs(jk - j0)) <= 1
        fprintf(1,'Terminated due to obj.func. changes \n');
        break;
    else
        jk = j0;
    end
    
    % check constraint violation
    if constraint >= tauStar
        break;
    end

end

return;                