% whether or not to show output
verbose = false;
verboseLevel = 0;

addpath ../optimsecond

nx = 45; ny = 45; nz = 1;
dx = 10; dy = 10; dz= 10;
G = cartGrid([nx ny nz],[nx*dx ny*dy nz*dz]);
G = computeGeometry(G);

rock.perm = txt2mat('two_frac2.out');
rock.poros = repmat(0.2, [G.cells.num, 1]);

fluid  = initTwoPhaseIncompressibleFluid('muo', 1, 'muw', 1,'swc',0.2,'sor',0.2);
resSolInit = initResSol(G, 0.0);
resSolInit.sw = .2*ones(G.cells.num, 1);

S = assembleMimeticSystem(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, 'InnerProduct', 'ip_tpf');
S.type = 'mixed';

% Choose objective function
objectiveFunction = str2func('NPV_brouwer');

% Hessian-vector product function
HessianFunction = str2func('HessMultFcn');

% Introduce wells
radius = .1;
W = [];
for i=1:45
   W = putWell(G, rock,W, (i-1)*nx +1      , 'Type', 'rate' ,'Val',  4  ,'Radius', radius, 'Name', ['I',num2str(i)]);
end

for i=1:45
   W = putWell(G, rock,W, i*nx      , 'Type', 'rate' ,'Val',  -4 ,'Radius', radius, 'Name', ['P',num2str(i)]);
end

W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Schedule: one time step, [0 1 PVI]
totVol   = sum(G.cells.volumes.*rock.poros);
schedule = initSchedule(W, 'TimeSteps', 120:120:240, 'Verbose', verbose);

injMaxMin = repmat( [0.001 200], 45, 1);
prodMaxMin = repmat( [-200 -0.001], 45, 1);


controls = initControls(schedule, 'ControllableWells', (1:90), ...
                                  'BHPMinMax', [0.001 200], ...
                                  'MinMax', [injMaxMin; prodMaxMin],...
                                  'Verbose', verbose);
                              
% % perform derivative checks
% usr_par = {G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction};
% Uinit = [controls.well.values]';
% Uinit = Uinit(:);
% deriv_check( Uinit, 1, usr_par);
                   
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

options = optimset('HessMult', HessianFunction, 'MaxIter', 30);

% Specify some model characteristics.
A   = []; 
b   = [];

Aeq = [ ones(1,90)   zeros(1,90);
        zeros(1,90)  ones(1,90) ];

beq = zeros(2,1);


tic;
% Call KNITRO to solve the optimization model.
[U,fval,exitflag,output,lambda] = ktrlink(@(U)NPVf_brouwer(U,auxdata),U,A,b,Aeq,beq,lb,ub,[],options, 'knitro.opt');

toc;