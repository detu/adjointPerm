%
%     function [ usr_par ] = xnew( u, iter, usr_par )
%
%     Purpose:
%
%     Compute the solution of the state equation and adjoint equation
%     for the optimal control of the unsteady Burger's equation. 
%     these solutions are stored in usr_par
%
%     Parameters
%  
%     On entry:
%
%     u      Variable  u. 
%
%     iter   Iteration counter if  new = 10  or  new = 20. 
%
%     usr_par user defined parameter. Used to pass problem
%            specific information.
%
%
%
%     On return:
%
%     usr_par user defined parameter. 
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


function [ usr_par ] = xnew( u, iter, usr_par )

% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.
global BURGERS_GLB

BURGERS_GLB.is_state_computed   = 0;
BURGERS_GLB.is_adjoint_computed = 0;

if( ~BURGERS_GLB.is_state_computed )
    % solve Burgers equation. Store the state as a global variable
    [ BURGERS_GLB.y, iflag ]      = state( u, BURGERS_GLB );
    BURGERS_GLB.is_state_computed = 1;
end
if( ~BURGERS_GLB.is_adjoint_computed )
     % solve the adjoint equation. Store the adjoint as a global variable
     [ BURGERS_GLB.lambda, iflag ]   = adjoint( BURGERS_GLB.y, BURGERS_GLB );
     BURGERS_GLB.is_adjoint_computed = 1;
end

%return


% plot the current control, the state and the adjoint
nt     = BURGERS_GLB.nt;       % number of time intervals
nx     = BURGERS_GLB.nx;       % number of spatial intervals
  
% plot the control
subplot(2,2,1)
mesh((0:BURGERS_GLB.Deltax:1), (0:BURGERS_GLB.Deltat:1), reshape(u,nx+1,nt+1)')
xlabel('x')
ylabel('t')
title(['Control u at iteration ',int2str(iter)])

% for plotting purposes augment y by the zero boundary condition
yplot = [ zeros(1,nt+1);
          BURGERS_GLB.y;
          zeros(1,nt+1)];
subplot(2,2,2)
mesh((0:BURGERS_GLB.Deltax:1), (0:BURGERS_GLB.Deltat:1), yplot')
xlabel('x')
ylabel('t')
title('Corresponding state y(u)')

% for plotting purposes augment lambda by the zero boundary condition
lplot = [ zeros(1,nt+1);
          BURGERS_GLB.lambda;
          zeros(1,nt+1)];
subplot(2,2,3)
mesh((0:BURGERS_GLB.Deltax:1), (0:BURGERS_GLB.Deltat:1), lplot')
xlabel('x')
ylabel('t')
title('Corresponding adjoint \lambda(u)')

drawnow;

%
% End of xnew.
