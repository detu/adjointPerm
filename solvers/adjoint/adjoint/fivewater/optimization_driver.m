% optimization_driver for Five water-spot case
%
% Synopsis
%   just type optimization_driver
%
% Description
%   Second order adjoint optimization using oil reservoir model
%   Newton conjugate method, will be explored other methods
%   such as: Lancos, Truncated Newton, preferable methods that utilize
%   distributed/parallel computing
% 
% Inputs ([]s are optional)
%   none
%  
%
% Outputs ([]s are optional)
%   none
%   
%
% Examples
%   type : optimization_driver
%
% See also
%   computeNumericalGradient (to check 1st gradient)
%
% Requirements
%   SINTEF Oil reservoir toolbox
%   can be downloaded from
%   http://www.sintef.no/Projectweb/GeoScale/Simulators/MRST/
%
% References
%   based on Matthias Heinkenschloss code
%   http://www.caam.rice.edu/~heinken/software/matlab_impl_constr/
%
% Authors
%   Eka Suwartadi (suwartad@itk.ntnu.no)
%
% License
%   The program is free for non-commercial academic use. Please 
%   contact the authors if you are interested in using the software
%   for commercial purposes. The software must not modified or
%   re-distributed without prior permission of the authors.
%
% Changes
%   03.05.2009  First development
%   ...      

clear all

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
addpath ../optimsecond
 
% initSimpleModel

% whether or not to show output
verbose = false;
verboseLevel = 0;

% load test/data/silvia.mat;

nx = 31; ny = 41; nz = 1;
cellDims = [nx ny nz];
physDims = [nx ny nz];

G = cartGrid(cellDims, physDims);

% set zero-poros cells inactive
% rock.poros = max(rock.poros, 1e-3);
rock.perm  = exp(( 3*rand(nx*ny, 1) + 1))*100  ;
rock.poros = repmat(0.3, [G.cells.num, 1]);

G = computeGeometry(G);

fluid  = initTwoPhaseIncompressibleFluid('muw', 0.3, 'muo', 3); %, 'swc', .2, 'sor', .2);
resSolInit    = initResSol(G, 0.0);
resSolInit.sw = .2*ones(G.cells.num, 1);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_simple');

% Choose objective function
objectiveFunction = str2func('simpleNPV');

% Hessian-vector product function
HessianFunction = str2func('HessMultFcn');

totRate = 2;
% Introduce wells
radius = .1;
W = putWell(G, rock,[], nx*20 + 16    , 'Type', 'rate' ,'Val',    -1 , 'Radius', radius, 'Name', 'p1');
W = putWell(G, rock, W, 1             , 'Type', 'bhp', 'Val',      1 , 'Radius', radius, 'Name', 'i1');
W = putWell(G, rock, W, nx            , 'Type', 'bhp', 'Val',      1 , 'Radius', radius, 'Name', 'i2');
W = putWell(G, rock, W, nx*ny-nx+1    , 'Type', 'bhp', 'Val',      1 , 'Radius', radius, 'Name', 'i3');
W = putWell(G, rock, W, nx*ny         , 'Type', 'bhp', 'Val',      1 , 'Radius', radius, 'Name', 'i4');

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

totVol   = sum(G.cells.volumes.*rock.poros);
schedule = initSchedule(W, 'NumSteps', 4, 'TotalTime', 0.5*totVol, 'Verbose', verbose);


controls = initControls(schedule, 'ControllableWells', (2:5), ...
                       'Verbose', verbose,...
                       'NumControlSteps', 4);


% perform derivative checks
usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
Uinit = [controls.well.values]';
Uinit = Uinit(:);
deriv_check( Uinit, 1, usr_par);


% Solve the optimal control problem using Newton's method
fprintf(1,'\n\n Solve the optimal control problem using')

% set initial iterate
u = Uinit;

% optimizer = 1;
% 
% if( optimizer == 1)
    fprintf(1,' Newton''s method\n')
    % set options for Newton CG method
    options.iprint = 1; % print level, default = 0
    options.maxit  = 100; %  maximum number of iterations, default = 100
    options.gtol   = 1.e-6; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-6; % stepsize stopping tolerance, default = 1.e-8
    
    tic;
    [ un, itern, iflagn ]  = newton_cg( u, options, usr_par);
    t = toc; fprintf(1, ' CPU time = %9.3f sec,   \n', t);
    fprintf(1,' Newton''s method returned after %d iterations with flag %d \n', itern,iflagn )

% elseif( optimizer == 2)
    fprintf(1,' limited memory BFGS method\n')
    % set options for limited memory BFGS method
    options.iprint = 1; % print level, default = 0
    options.gtol   = 1.e-6; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-6; % stepsize stopping tolerance, default = 1.e-8
    options.maxit  = 100; % Maximum number of iterations, default 100
    options.L      = 20;  % LBFGS storage limit, default 20
    tic;
    [ ul, iterl, iflagl ]  = lbfgs( u, options, usr_par);
    t = toc; fprintf(1,' CPU time = %9.3f sec,   \n', t);
    fprintf(1,' LBFGS returned after %d iterations with flag %d \n', iterl,iflagl )
    
%     fprintf(1,' limited memory BFGS method\n')
%     % set options for limited memory BFGS method
%     options.iprint = 1; % print level, default = 0
%     options.gtol   = 1.e-5; % gradient stopping tolerance, default = 1.e-8
%     options.stol   = 1.e-5; % stepsize stopping tolerance, default = 1.e-8
%     options.maxit  = 20; % Maximum number of iterations, default 100
%     options.L      = 20;  % LBFGS storage limit, default 20
%     tic;
%     [ ul1, iterl1, iflagl1 ]  = lbfgs1( u, options, usr_par);
%     t = toc; fprintf(1,' CPU time = %9.3f sec,   \n', t);
%     fprintf(1,' LBFGS returned after %d iterations with flag %d \n', iterl1,iflagl1 )

% else
%     fprintf(1,' optimizer is %d; must be 1, 2 \n', optimizer)
% end





