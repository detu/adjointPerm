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
% addpath ../optimization
 
% initSimpleModel

% whether or not to show output
verbose = false;
verboseLevel = 0;
nx = 5; ny = 5; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
m     = rand(5,5);

Kreal = [ 1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1; ...
          1 1 1 1 1 ];
      
m     = reshape(m, [1,G.cells.num]);
Kreal = reshape(Kreal, [1,G.cells.num]);
K     = Kreal ./ m;

% save Kreal.mat Kreal m;
save Kreal.mat Kreal;
perm      = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

% Choose objective function
objectiveFunction = str2func('FwMatchSmallReg');
% objectiveFunction = str2func('FwMatchSmall');

modelHMSmall3

% totVol   = sum(G.cells.volumes.*rock.poro);
schedule = initSchedule(W, 'NumSteps', 1, 'TotalTime', 100*day, 'Verbose', verbose);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', 1);


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





