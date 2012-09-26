% whether or not to show output
verbose = false;
verboseLevel = 0;

load test/data/silvia.mat;
addpath ../optimsecond

nx = 31; ny = 41; nz = 1;
cellDims = [nx ny nz];
physDims = [nx ny nz]*2;

G = cartGrid(cellDims, physDims);

% set zero-poros cells inactive
rock.poros = repmat(0.2, [G.cells.num, 1]);
rock.perm  = 1e-2*rock.perm;

G = computeGeometry(G);

fluid  = initTwoPhaseIncompressibleFluid('muw', 0.3, 'muo', 3, 'swc', .2, 'sor', .2);
resSolInit    = initResSol(G, 0.0);
resSolInit.sw = .2*ones(G.cells.num, 1); 

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');

% Choose objective function
objectiveFunction = str2func('simpleNPV');

% Hessian-vector product function
HessianFunction = str2func('HessMultFcn');

totRate = 10;
% totRate = 40;
% Introduce wells
radius = .1;
W = putWell(G, rock,[], nx*20 + 16    , 'Type', 'rate' ,'Val',     totRate , 'Radius', radius, 'Name', 'p1');
W = putWell(G, rock, W, 1             , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i1');
W = putWell(G, rock, W, nx            , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i2');
W = putWell(G, rock, W, nx*ny-nx+1    , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i3');
W = putWell(G, rock, W, nx*ny         , 'Type', 'rate', 'Val',-.25*totRate , 'Radius', radius, 'Name', 'i4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poros);
schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', 0.2*totVol, 'Verbose', verbose);


injMaxMin = repmat( [0.001 100], 1, 1);
prodMaxMin = repmat( [-100 -0.001], 4, 1);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'LinEqConst', {ones(1,5)}, ..., 
                       'Verbose', verbose, ...
                       'NumControlSteps', 5);

% perform derivative checks
usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
Uinit = [controls.well.values]';
Uinit = Uinit(:);
deriv_check( Uinit, 1, usr_par);
                   
% % -------Setting of IPOPT inputs------------
% % initial  control input
% U = [controls.well.values]';
% U = U(:);
% 
% minMax = vertcat(controls.well.minMax);
% numstep = numel(controls.well);
% 
% minMax = repmat(minMax,numstep,1);
% lb = minMax(:,1);
% ub = minMax(:,2);
% 
% lbc = [0 0 0 0 0];
% ubc = [0 0 0 0 0];
% 
% % adjusting barrier parameter to objective function
% auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
% HessianFunction = @(varargin) HessMultFcn(varargin{:}, auxdata);
% 
% options = optimset('HessMult', HessianFunction, 'MaxIter', 30);
% 
% A   = []; 
% b   = [];
% Aeq = [ones(1,5) zeros(1,20);... 
%        zeros(1,5) ones(1,5) zeros(1,15);...
%        zeros(1,10) ones(1,5) zeros(1,10);...
%        zeros(1,15) ones(1,5) zeros(1,5);...
%        zeros(1,20) ones(1,5) ];
%    
% beq = zeros(5,1);
%              
% tic;
% % run KNITRO
% [U,fval,exitflag,output,lambda] = ktrlink(@(U)NPV_hess(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options, 'knitro.opt');
% toc;

% Solve the optimal control problem using Newton's method
fprintf(1,'\n\n Solve the optimal control problem using')

% set initial iterate
u = Uinit;

% 
% if( optimizer == 1)
    fprintf(1,' Newton''s method\n')
    % set options for Newton CG method
    options.iprint = 1; % print level, default = 0
    options.maxit  = 20; %  maximum number of iterations, default = 100
    options.gtol   = 1.e-2; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-2; % stepsize stopping tolerance, default = 1.e-8
    
    tic;
    [ un, itern, iflagn ]  = newton_cg( u, options, usr_par);
    t = toc; fprintf(1, ' CPU time = %9.3f sec,   \n', t);
    fprintf(1,' Newton''s method returned after %d iterations with flag %d \n', itern,iflagn )

% elseif( optimizer == 2)
    fprintf(1,' limited memory BFGS method\n')
    % set options for limited memory BFGS method
    options.iprint = 1; % print level, default = 0
    options.gtol   = 1.e-2; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-2; % stepsize stopping tolerance, default = 1.e-8
    options.maxit  = 20; % Maximum number of iterations, default 100
    options.L      = 20;  % LBFGS storage limit, default 20
    tic;
    [ ul, iterl, iflagl ]  = lbfgs( u, options, usr_par);
    t = toc; fprintf(1,' CPU time = %9.3f sec,   \n', t);
    fprintf(1,' LBFGS returned after %d iterations with flag %d \n', iterl,iflagl )
    

