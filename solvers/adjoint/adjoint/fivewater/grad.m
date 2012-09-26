%     function [ fg ] = grad( u, usr_par )
%
%     Purpose:
%
%     Compute the gradient of the objective function f(y(u),u)
%     for the optimal control of the unsteady Burger's equation. 
%     In this implementation we assume equidistant time steps 
%     and g = 0.
%
%
%     Parameters
%
%     On entry:
%
%     u      Variable  u.
%            u((i-1)*(nx+1)+1), ..., u(i*(nx+1))
%            controls at time (i-1)*Deltat, i = 1, ..., nt+1
%
%     usr_par user defined parameter. Used to pass problem
%            specific information.
% 
%     On return:
%
%     fg      Value of the gradient of the objective function  f(y(u),u).
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

function [ fg ] = grad( u, usr_par )

% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.

global BURGERS_GLB

% get problem data
Deltat = BURGERS_GLB.Deltat;   % length of time interval
dt2    = Deltat/2;
Deltax = BURGERS_GLB.Deltax;   % length of spatial interval
nt     = BURGERS_GLB.nt;       % number of time intervals
nx     = BURGERS_GLB.nx;       % number of spatial intervals
ny     = nx-1;             % number of y-variables per time step
nu     = nx+1;             % number of u-variables per time step

omega  = BURGERS_GLB.omega;

if( ~BURGERS_GLB.is_state_computed )
    % solve Burgers equation. Store the state as a global variable
    [ BURGERS_GLB.y, iflag ]      = state( u, BURGERS_GLB );
    BURGERS_GLB.is_state_computed = 1;
end
if( ~BURGERS_GLB.is_adjoint_computed )
     % solve the adjoint equation. Store the adjoint as a global variable
     [ BURGERS_GLB.lambda, iflag ]    = adjoint( BURGERS_GLB.y, BURGERS_GLB );
     BURGERS_GLB.is_adjoint_computed  = 1;
end

  
% allocate space for the gradient
fg = zeros(size(u));
  
  
fg(1:nu) = (dt2*omega)*(BURGERS_GLB.Q*u(1:nu)) ...
                  + dt2*(BURGERS_GLB.B'*BURGERS_GLB.lambda(:,2));

for i = 2:nt
	fg((i-1)*nu+1:i*nu) = (Deltat*omega)*(BURGERS_GLB.Q*u((i-1)*nu+1:i*nu))...
                               + dt2*(BURGERS_GLB.B'*(BURGERS_GLB.lambda(:,i)...
                                                  +BURGERS_GLB.lambda(:,i+1)));
end
i = nt+1;
fg((i-1)*nu+1:i*nu) = (dt2*omega)*(BURGERS_GLB.Q*u((i-1)*nu+1:i*nu))...
                                     + dt2*(BURGERS_GLB.B'*BURGERS_GLB.lambda(:,i));
  
% End of grad.



