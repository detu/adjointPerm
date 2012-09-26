% optimization_driver for Five-spot case
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

% whether or not to show output
verbose = false;
verboseLevel = 0;

[G, S, rock, W, fluid] = initModel(2, 1);
% [G, S, rock, W, fluid] = initModel1(1, 1);

resSolInit    = initResSol(G, 0.0);
resSolInit.sw = 0.2*ones(G.cells.num, 1);


% Choose objective function
% objectiveFunction = str2func('simpleNPV');
objectiveFunction = str2func('simpleNPVP_tid');

% Hessian-vector product function
HessianFunction = str2func('HessMultFcn');

% totRate = 2;

totVol   = sum(G.cells.volumes.*rock.poros);
schedule = initSchedule(W, 'NumSteps', 5, 'TotalTime', 2*totVol, 'Verbose', verbose);

injMaxMin  = repmat( [0.001 1000000], 1);
prodMaxMin = repmat( [-1000000 -0.001], 4, 1);

controls = initControls(schedule, 'ControllableWells', (1:5), ...
                       'MinMax', [injMaxMin; prodMaxMin],...
                       'Verbose', verbose,...
                       'NumControlSteps', 5);

% controls = initControls(schedule, 'ControllableWells', (1:5), ...
%                        'Verbose', verbose,...
%                        'NumControlSteps', 5);


% % perform derivative checks
usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
Uinit = [controls.well.values]';
Uinit = Uinit(:);
% deriv_check( Uinit, 1, usr_par);
% return;

%---OPTIMIZATION PART---------
% initial  control input
U = [controls.well.values]';
U = U(:);
minMax = vertcat(controls.well.minMax);
numstep = numel(schedule);
minMax = repmat(minMax,numstep,1);
lb = minMax(:,1);
ub = minMax(:,2);

auxdata = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
HessianFunction = @(varargin) HessMultFcn(varargin{:}, auxdata);

options = optimset('HessMult', HessianFunction, 'MaxIter', 50);

% Specify some model characteristics.
A   = []; 
b   = [];

Aeq = [ ones(1,5)   zeros(1,20);
        zeros(1,5)  ones(1,5)  zeros(1,15);
        zeros(1,10) ones(1,5)  zeros(1,10);
        zeros(1,15) ones(1,5)  zeros(1,5);
        zeros(1,20) ones(1,5) ];

beq = zeros(5,1);
% Aeq = [];
% beq = [];


tic;
% Call KNITRO to solve the optimization model.
[U,fval,exitflag,output,lambda] = ktrlink(@(U)NPVf_brouwer(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options, 'knitro.opt');

toc;

% Ubfgs = U;
% save BFGSResults.mat;

Utn = U;
save TNResults.mat;