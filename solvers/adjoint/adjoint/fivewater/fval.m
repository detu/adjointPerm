%     function [ f ] = fval( u, usr_par )
%
%     Purpose:
%
%     Compute the value of the objective function f(y(u),u)
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
%     f      Value of the objective function  f(y(u),u).
%
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

  function [ f ] = fval( u, usr_par )
  
  % We use global variables to pass information to the
  % functions that compute the solution of Burgers equation, etc.
  global BURGERS_GLB
  
  % get problem data
  Deltat = BURGERS_GLB.Deltat;   % length of time interval
  Deltax = BURGERS_GLB.Deltax;   % length of spatial interval
  nt     = BURGERS_GLB.nt;       % number of time intervals
  nx     = BURGERS_GLB.nx;       % number of spatial intervals
  nx1    = nx+1;
  omega  = BURGERS_GLB.omega;    % control penalty

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

  
  f = (Deltat/4)*(BURGERS_GLB.y(:,1)'*BURGERS_GLB.M*BURGERS_GLB.y(:,1)) ...
       + (Deltat/2)*(BURGERS_GLB.y(:,1)'*BURGERS_GLB.g(:,1)) ...
       + (Deltat*omega/4)*(u(1:nx1)'*BURGERS_GLB.Q*u(1:nx1));
  for i = 2:nt
      f = f + (Deltat/2)*(BURGERS_GLB.y(:,i)'*BURGERS_GLB.M*BURGERS_GLB.y(:,i)) ...
            + Deltat*(BURGERS_GLB.y(:,i)'*BURGERS_GLB.g(:,i)) ...
            + (Deltat*omega/2)*(u((i-1)*nx1+1:i*nx1)'*BURGERS_GLB.Q*u((i-1)*nx1+1:i*nx1));
  end
  i = nt+1;
  f = f + (Deltat/4)*(BURGERS_GLB.y(:,i)'*BURGERS_GLB.M*BURGERS_GLB.y(:,i)) ...
        + (Deltat/2)*(BURGERS_GLB.y(:,i)'*BURGERS_GLB.g(:,i)) ...
        + (Deltat*omega/4)*(u((i-1)*nx1+1:i*nx1)'*BURGERS_GLB.Q*u((i-1)*nx1+1:i*nx1));

  
% End of fval.



