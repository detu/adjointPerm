function [outputArgs] = setupCaseHM1( )
% History Matching Case
% Implementation of Hessian-Free Newton method 
% Aim: 1. To compare quasi Newton methods: BFGS, SR1 against Newton-CG
%      2. To assess the advantage of Line search methods vs. Trust Region
%      based methods 
% 
% [OUTPUTARGS] = SETUPCASEHM1(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2010/05/27 15:16:56 $	$Revision: 0.1 $
% Copyright: IO Center NTNU - Norway 2010

% reservoir properties : 
% number of grid block dimension
nx       = 49;
ny       = 49;
nz       = 1;

% number of parameter, 
% in this case permeability, 
% which will be matched
numParam = nx * ny * nz; 

% initial parameter values 
% for optimization, i.e. 1 mD(unit)
x0       = ones(numParam,1);

% objective function calculation
% and gradient computation
fn       = @mismatchHM;

% % LOAD STARTFILE.M FIRST ???
% startfile

opts = lbfgs_options('iprint', 1, 'maxits', 1000, 'factr', 1e5, ...
                     'cb', @test_callback);
                 
[x,fx,exitflag,userdata] = lbfgs(fn,x0,opts);

% SAVE PERM ???

semilogy(userdata.f);
xlabel('Iteration');
ylabel('f');

end
