function [] = prob_gen(nx, nt, viscosity, omega)
%  
% Set problem data for the optimal control of Burgers' equation
%
% inputs:
%    nx         number of spatial intervals used in the discretization
%               The spatial mesh points are   
%                      0, 1/nx, 2/nx, ..., (nx-1)/nx, 1
%
%    nt         number of temporal intervals used in the discretization
%               The time steps are
%                     0, 1/nt, 2/nt, ..., (nt-1)/nt, 1
%
%    viscosity  viscosity parameter
%
%    omega      penalty parameter for the control
%
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%
% 

% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.

global BURGERS_GLB

% set problem data
BURGERS_GLB.nt         = nt;
BURGERS_GLB.nx         = nx;
BURGERS_GLB.Deltat     = 1/BURGERS_GLB.nt;   % length of time interval
BURGERS_GLB.Deltax     = 1/BURGERS_GLB.nx;   % length of spatial interval
BURGERS_GLB.viscosity  = viscosity;             % viscosity
BURGERS_GLB.omega      = omega;          % control penalty parameter

Deltax = 1/nx;  

e = ones(nx-1,1);
% mass matrix
BURGERS_GLB.M = (Deltax/6)*spdiags([e 4*e e], -1:1, nx-1, nx-1);
% stiffness matrix
BURGERS_GLB.A = (viscosity/Deltax)*spdiags([-e 2*e -e], -1:1, nx-1, nx-1);

e = ones(nx+1,1);
B = -(Deltax/6)*spdiags([e 4*e e], -1:1, nx+1, nx+1);
BURGERS_GLB.B = B(2:nx,:);

% 
BURGERS_GLB.Q            = (Deltax/6)*spdiags([e 4*e e], -1:1, nx+1, nx+1);
BURGERS_GLB.Q(1,1)       = Deltax/3;
BURGERS_GLB.Q(nx+1,nx+1) = Deltax/3;

% vectors r in the righ hand side of Burgers equation
% r(:,i)  vector r at time step i
BURGERS_GLB.r          = zeros(nx-1,nt+1);


% vectors g in the objective function
% g(:,i)  vector g at time step i
BURGERS_GLB.g          = zeros(nx-1,nt+1);
nx2                    = floor(nx/2);
BURGERS_GLB.g(1:nx2,:) = -BURGERS_GLB.Deltax*ones(nx2,nt+1);
