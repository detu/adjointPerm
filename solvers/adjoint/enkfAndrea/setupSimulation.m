function [Xt0, Xen0, U0, auxdata, R, H] = setupSimulation(Ne,trueRealization)

% Inputs
% Ne : is the number of realizations in the ensemble
% trueRealization: is the number of true realization in the ensemble

% Outputs:
% Xt0 : is the initial augmented state for the true realization, X = [fw; K; Sw; P]
% Xen0 : is the initial augmented state for the ensemble realizations
% U0 : is the initial input vector
% R: observation error covariance
% H: observation operator, Yt = H X
% optimal control U

%% CASE DESCRIPTION (to update)
% A simple example of well-rate control
% 5-spot pattern: 4 producer + 1 injector 
% reservoir geometri taken from layer 36th SPE 10th
% 2D: 31 x 41 gridblocks
% physical dimension of a grid block 100 ft \times 100 ft \times 100 ft
% mobility ratio oil-water: 3:0.3
% end point: Sor, Swc: 0.2 : 0.2
% initial water saturation: 0.2
% objective function: NPV (see the function def. for detailed prices)
% initial control: 5000 ft^3/day injector

%% PATH INITIALIZATION

% % 
% % %cd /home/acap/programs/EKA_MRST;
% % cd C:\Users\acap\Documents\phdRepository\Trondheim\EKA_MRST
% % 
% % startup;
% % %cd /home/acap/programs/EKA_MRST/solvers/adjoint/tests/myTests;
% % cd C:\Users\acap\Documents\phdRepository\Trondheim\EKA_MRST\solvers\adjoint\tests\myTests
% % 
opengl neverselect;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);

% whether or not to show output


%% RESERVOIR SETUP
%load realizations.mat; 
%load myPermRob.mat;

% load spe realizations
load myPermRob.mat;


nx = 45%31; 
ny = 45%41; 
nz = 1;
cellDims = [nx ny nz];
physDims = [10*nx 10*ny 10*nz];



%% Reservoir parameters

G = cartGrid(cellDims, physDims);
G = computeGeometry(G);





% model parametes
modelPar.verbose = false;
modelPar.verboseLevel = 0;

rock.poro  = 0.2;
fluid      = initCoreyFluid('mu', [0.3 3], 'sr', [.1 .1]);

wellPar.injIdx = nx*ceil(ny/2) + ceil(nx/2); %indexes of iniectors
wellPar.prodIdx = [1, (nx), (nx*ny-nx+1), (nx*ny)] ; %indexes of producers

wellPar.wellRadius = 0.1; % meters

% simulation parameters
simParam.simTime = 810;%1090;% time horizon in days
simParam.numContrSteps = 6;%; number of control steps in the time horizon
simParam.totRate = 500/day; % cubic meters a day
simParam.objectiveFunction = str2func('my_NPV'); % Choose objective function
%objectiveFunction = str2func('recovery');


poreVol = mean(rock.poro)*prod(physDims); % pore volume
numDayPore = poreVol/(simParam.totRate* day);
disp(['numbers of days required to inject one pore volume: ' num2str(numDayPore)])




auxdata = {G, rock, fluid,wellPar,simParam,modelPar};


%% generate initial augmented state vector X = [fw; K; Sw; P] for both true and ensemble states
totRate = simParam.totRate;
prodIdx = wellPar.prodIdx;
injIdx = wellPar.injIdx;

numInj = length(injIdx);
numProd = length(prodIdx);
numContrSteps = simParam.numContrSteps;

Xt0 = zeros( (numProd + 3*prod(cellDims) ), 1 );
Xen0 = zeros( (numProd + 3*prod(cellDims) ), Ne );

fw = zeros(length(prodIdx),1);
perm = pReal(:,trueRealization)*darcy;
P = 300*barsa() * ones(prod(cellDims),1); % initial pressure in bar (300 bar)
S = 0.1 * ones(prod(cellDims),1); % initial oil saturation


Xt0 = [fw; perm; S; P];

% ensemble realizations
for k1=1:Ne
    perm = pReal(:,k1)*darcy;
    Xen0(:,k1) = [fw; perm; S; P];
end

%% generate initial input (equals for true and ensemble realizations)
U0 = repmat( (totRate/numInj), numInj,1 );
U0 = [U0; repmat( -(totRate/numProd), numProd ,1 )];

U0 = repmat(U0,numContrSteps,1);



%% observation error covariance
R = 1e-2;
%% observation operator 
H = [speye(numProd), sparse(numProd, (length(Xt0)- numProd))];     

