%
%     function [ y, iflag ] = state( u, usr_par )
%
%     Purpose:
%
%     Given u, solve the unsteady Burgers equation c(y,u) = 0.
%
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
%    
%     On return:
%
%     y      solution of  c(y,u) = 0.
%            y(:,i)  solution at time  (i-1)*Deltat, i = 1, ..., nt+1
%
%     iflag  Error flag
%            iflag = 0 if the computation was successful;
%            iflag > 0 otherwise.
%
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

  function [ y, iflag ] = state( u, usr_par )
  
  % We use global variables to pass information to the
  % functions that compute the solution of Burgers equation, etc.
  global BURGERS_GLB
  
  iflag = 0;
  idebg = 0; % set idebg > 0 to get convergence history of Newton's method
  
  % get problem data
  Deltat = BURGERS_GLB.Deltat;   % length of time interval
  Deltax = BURGERS_GLB.Deltax;   % length of spatial interval
  nt     = BURGERS_GLB.nt;       % number of time intervals
  nx     = BURGERS_GLB.nx;       % number of spatial intervals
  ny     = nx-1;             % number of y-variables per time step
  nu     = nx+1;             % number of u-variables per time step

  dt2    = Deltat/2;
  
  % initial time
  rhs        = zeros(nx-1,1);
  nx2        = floor(nx/2);
  rhs(1:nx2) = ones(nx2,1);
  y(:,1)     = rhs;
  %y(:,1)     = BURGERS_GLB.M\(Deltax*rhs);
  
  for i = 1: nt
      % time step i
      
      % compute y(:,i+1) using Newton's method
      tol      = 1.e-2*min(Deltat^2,Deltax^2);
      max_iter = 10;

      % initial guess 
      y(:,i+1) = y(:,i);
      
      % Residual (res0 is the residual component that is independent of y(:,i+1))
      res0 = (-BURGERS_GLB.M+dt2*BURGERS_GLB.A)*y(:,i) + dt2*Ny(y(:,i)) ...
             + dt2*(BURGERS_GLB.B*(u((i-1)*nu+1:i*nu) + u(i*nu+1:(i+1)*nu))) ...
             - dt2*(BURGERS_GLB.r(:,i) + BURGERS_GLB.r(:,i+1));
      
      res  = res0 + (BURGERS_GLB.M+dt2*BURGERS_GLB.A)*y(:,i+1) + dt2*Ny(y(:,i+1)); 
                     
      iter = 0;
      
      if( idebg > 0 )
         disp('   time step    Newton iter   Residual    Step size'); format short e
      end
      
      while (iter < max_iter & norm(res) >= tol )
      
           resnorm = norm(res);
           if( idebg > 0 )
             if( iter == 0 )
                 disp([ i    iter     resnorm ])
             else
                 disp([ i    iter     resnorm  stepsize ])
             end
           end
           
           
           % Jacobian
           Mat = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(y(:,i+1)));
           
           % Newton step
           stepy = - Mat\res;
           
           % Compue new guess using Armijo rule
           stepsize = 1;         % step size
           ytmp     = y(:,i+1) + stepsize*stepy;
           res      = res0 + (BURGERS_GLB.M+dt2*BURGERS_GLB.A)*ytmp + dt2*Ny(ytmp);
           
           while ( norm(res) >= (1-1.e-4*stepsize)*resnorm )
               stepsize = stepsize/2;
               ytmp     = y(:,i+1) + stepsize*stepy;
               res      = res0 + (BURGERS_GLB.M+dt2*BURGERS_GLB.A)*ytmp + dt2*Ny(ytmp);
           end
           y(:,i+1) = ytmp;
           
           iter = iter+1;
      end
  end
           
%
% End of state.
