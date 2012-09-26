%
%     function [ lambda, iflag ] = adjoint( y, usr_par )
%
%     Purpose:
%
%     Compute  lambda  by solving the adjoint equation for the
%     unsteady Burger's equation. In this implementation we
%     assume equidistant time steps and g = 0.
%
%   
%     Parameters
%
%     On entry:
%
%     y      Variable  y.
%            y(:,i)  solution of Burgers equation
%            at time  (i-1)*Deltat, i = 1, ..., nt+1
%
%     usr_par user defined parameter. Used to pass problem
%            specific information.
%
%     On return:
%
%     lambda Solution of the adjoint equation.
%            lambda(:,i)  solution of adjoint equation
%            at time  (i-1)*Deltat, i = 1, ..., nt+1
%
%
%     iflag  Error flag.
%            iflag = 0 if the computation was successful;
%            iflag > 0 otherwise.
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%
%
  function [ lambda, iflag ] = adjoint( y, usr_par )
  
  % We use global variables to pass information to the
  % functions that compute the solution of Burgers equation, etc.
  global BURGERS_GLB

  iflag = 0;
  
  % get problem data
  Deltat = BURGERS_GLB.Deltat;   % length of time interval
  dt2    = Deltat/2;
  Deltax = BURGERS_GLB.Deltax;   % length of spatial interval
  nt     = BURGERS_GLB.nt;       % number of time intervals
  nx     = BURGERS_GLB.nx;       % number of spatial intervals
  
  ny     = nx-1;  % number of y-variables per time step
  nu     = nx+1;   % number of u-variables per time step


  % allocate space for the adjoint lambda
  lambda  = zeros(ny,nt+1); 
  
  
  % solve for lambda
  Mat            = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(y(:,nt+1)));
  rhs            = -dt2*(BURGERS_GLB.M*y(:,nt+1) + BURGERS_GLB.g(:,nt+1));
  lambda(:,nt+1) = Mat'\rhs; 
    
  for i = nt:-1:2
     Mat            = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(y(:,i)));
     rhs            = -(-BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(y(:,i))))'*lambda(:,i+1)...
                       - Deltat*(BURGERS_GLB.M*y(:,i) + BURGERS_GLB.g(:,i));
     lambda(:,i)    = Mat'\rhs;
  end
%      
% End of adjoint.










