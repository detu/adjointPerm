function [Hv] = H0vec(v, usr_par)
% 
% compute H0*v, where H0 is the initial BFGS matrix
% H0 is a replacement of the inverse of the Hessian of f!
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.

% global BURGERS_GLB


% omega  = BURGERS_GLB.omega;
% u = 0;
% Hv   = -Hessvec( v, u, usr_par );
% Hv = zeros(size(v));
% Hv(1,1) = v(1,1); 
% Hv = v;
Hv     = ones(size(v));
return;

% use omega*I  as initial Hessian approx.
%Hv = v/omega;
%return;
    
% % use omega*dt*Q  as initial Hessian approx.
% Deltat = BURGERS_GLB.Deltat;   % length of time interval
% dt2    = Deltat/2;
% nt     = BURGERS_GLB.nt;       % number of time intervals
% nx     = BURGERS_GLB.nx;       % number of spatial intervals
% nu     = nx+1;                 % number of u-variables per time step
% 
% Hv(1:nu) = (omega*dt2)\(BURGERS_GLB.Q\v(1:nu));
% for i = 2:nt
%         Hv((i-1)*nu+1:i*nu) = (omega*Deltat)\(BURGERS_GLB.Q\v((i-1)*nu+1:i*nu));
% end
% i = nt+1;
% Hv((i-1)*nu+1:i*nu) = (omega*dt2)\(BURGERS_GLB.Q\v((i-1)*nu+1:i*nu));
    
