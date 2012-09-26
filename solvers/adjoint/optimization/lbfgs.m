%
% function [ x, iter, iflag ]  = lbfgs( x, options, usr_par)
%
% Limited memory BFGS method for the solution of    min f(x)   
%
% The user has provide the functions:
%     function [fx] = fval(x, usr_par)   evaluate f(x)
%     function [gx] = grad(x, usr_par)   evaluate grad_f(x)
%
%     function [Hv] = H0vec(v, usr_par)  compute H0*v, where
%                                        H0 is the initial BFGS matrix
%                                        (H0 ~ inverse of Hessian of f at
%                                        intial x).
%
%     function [usr_par] = xnew( x, iter, usr_par)
%     function [x1x2]    = xprod( x1, x2, usr_par).
%     
% xnew is called whenever a new x is generated and it is called
% before any of the three functions fval grad, Hessvec are called 
% with this new x.
%
% xprod evaluates the inner product of x1 and x2.
%
%
% input
%   x       real(n) 
%           initial iterate
%
%   options   Newton-CG parameters
%             options.iprint   print if ~=0, default = 0
%             options.fid      file identifier, default = 1
%             options.maxit    maximim number of iterations, default = 100
%             options.gtol     gradient stopping tolerance, default = 1.e-8
%             options.stol     stepsize stopping tolerance, default = 1.e-8
%             options.L        number of updates stored, default = 20
%             if lbfgs  is called with options =[], then defauls
%                   are used.
%         
%   usr_par problem specific information. usr_par is not referenced in
%           lbfgs, but passed to the user provided fucntions fval, grad, ....
%
%
% output:
%    x      final approximate solution
%    iter   number of Newton iterations
%    iflag  error flag
%           iflag = 0, Successful  norm(g) < gtol
%           iflag = 1, linesearch returned with error
%           iflag = 2, reached max number of iterations without achieving
%                      gradient norm > gtol
%           iflag = 3, step norm < stol, but gradient norm > gtol
%                 
%    
%   Matthias Heinkenschloss
%   Department of Computational and Applied Mathematics
%   Rice University
%   June 6, 2008
%
%
function [ x, iter, iflag ]  = lbfgs( x, options, usr_par)


%% set linear search parameters
iprint  = 0;
fid     = 1;
maxit   = 100;
gtol    = 1.e-8;
stol    = 1.e-8;
L       = 20;
if( ~isempty(options) )
    if( isfield(options,'iprint') ); iprint = options.iprint; end
    if( isfield(options,'fid') );    fid    = options.fid; end
    if( isfield(options,'maxit') );  maxit  = options.maxit; end
    if( isfield(options,'gtol') );   gtol   = options.gtol; end
    if( isfield(options,'stol') );   stol   = options.stol; end
    if( isfield(options,'L') );      L      = options.L; end
end

%% initialization
iflag       = 0;
iter        = 0;

% storage for vectors s, y, and scalars rho that define BFGS matrices
BFGSs   = zeros(length(x),L);
BFGSy   = zeros(length(x),L);
BFGSrho = zeros(L,1);
BFGSa   = zeros(L,1);
Lc      = 0;                  % current number of updates stored.



if( iprint > 0 )
   fprintf(1,' LBFGS \n')
   fprintf(1,['     k    f(x_k)    ||gradf(x_k)||    ||s_k||      stepsize\n'])
end

usr_par  = xnew( x, iter, usr_par);
f        = fval(x, usr_par);
g        = grad(x, usr_par);
gnorm    = sqrt(xprod( g, g, usr_par));
stepsize = 1;
snorm    = 2*stol;

%stopping criteria for BFGS
while ( ~(gnorm < gtol | stepsize*snorm < stol ) ...
        & iter <= maxit & iflag == 0)
   
   % Compute search direction  s = - H*g
   q = -g;
   for i = Lc:-1:1
       BFGSa(i) = BFGSrho(i) * xprod( BFGSs(:,i), q, usr_par);
       q        = q - BFGSa(i)*BFGSy(:,i);
   end
   s = H0vec(q, usr_par);
   for i = 1:Lc
       scal = xprod( BFGSy(:,i), s, usr_par);
       s    = s + (BFGSa(i) - BFGSrho(i)*scal)*BFGSs(:,i);
   end
   
   % check descent property
   gs = xprod( g, s, usr_par );
   
   if( gs >= 0 )
       if( iprint > 0 )
          fprintf(fid,' BFGS iteration %d \n', iter)
          fprintf(fid,' Direction s is not a descent direction')
          fprintf(fid,'  <grad,s> = %12.6e; terminate \n', gs)
       end
       iflag = 2;
       return
   end
   snorm = sqrt(xprod( s, s, usr_par));

   % compute step size
   if( iprint <= 1 )
       ln_options.iprint = 0;
   else
       ln_options.iprint = 1;
   end
   ln_options.stpmin = 0.1*stol/snorm;
   
   % save current iterate, fucntion value, gradient
   f_c      = f;  
   g_c      = g;  
   gnorm_c  = gnorm;
   x_c      = x;  
   stepsize = 1;  % initial stepsize
   
   [stepsize, x, f, iter_ls, iflag_ls, usr_par] ...
            = lnsrch_bt( x, s, gs, f_c, stepsize, ln_options, usr_par );
   g        = grad(x, usr_par);
       
   if( iflag_ls > 0 )
      if( iprint > 0 )
          fprintf(fid,' LBFGS iteration %d \n', iter)
          fprintf(fid,' line search returned with flag = %d; terminate \n', iflag_ls)
      end
      iflag = 1;
      return
   end
   
   gnorm = sqrt(xprod( g, g, usr_par));
   
   if( iprint > 1 )
   fprintf(fid,' lbfgs \n')
   fprintf(fid,['     k    f(x_k)    ||gradf(x_k)||    ||s_k||       stepsize\n'])
   end
   if( iprint > 0 )
      fprintf(fid,'  %4d  %12.6e  %12.6e  %12.6e  %12.6e \n', ...
                  iter, f_c, gnorm_c, snorm, stepsize )
   end
   
   % update BFGS marices
   if ( Lc == L )
       % storage limit reached
       BFGSs(:,1:L-1) = BFGSs(:,2:L);
       BFGSy(:,1:L-1) = BFGSy(:,2:L);
       BFGSrho(1:L-1) = BFGSrho(2:L);
       Lc             = L-1;
   end
   Lc = Lc+1;
   BFGSs(:,Lc) = x - x_c;
   BFGSy(:,Lc) = g - g_c;
   BFGSrho(Lc) = 1/xprod( BFGSy(:,Lc), BFGSs(:,Lc), usr_par);
   if( BFGSrho(Lc) <= 0 )
      if( iprint > 0 )
          fprintf(fid,' LBFGS iteration %d \n', iter)
          fprintf(fid,' rho  = %12.5e < 0; skip update \n', BFGSrho(Lc))
      end
      Lc = Lc-1;
   end
   
       

   iter  = iter+1;
end

if( iprint > 0 & iflag == 0 )
      fprintf(fid,'  %4d  %12.6e  %12.6e \n\n', iter, f, gnorm )
end

if ( gnorm > gtol & iter >= maxit)
   iflag = 2;
end
if ( gnorm > gtol & stepsize*snorm >= stol)
   iflag = 3;
end

   
   
   



