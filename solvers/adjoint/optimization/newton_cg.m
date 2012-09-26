%
% function [ x, iter, iflag ]  = newton_cg( x, options, usr_par)
%
% Newton Conjugate Gradient method for the solution of    min f(x)   
%
% newton_cg computes the search direction using the CG method applied
% to    min < gradf(x), s > + 0.5 * < hessf(x)*s, s >
% and it computes the new iterate  x + stepsize* s,
% where the step size is computed using a backtracking linear search 
% procedure.
%
% The user has provide the functions:
%     function [fx] = fval(x, usr_par)         evaluate f(x)
%     function [gx] = grad(x, usr_par)         evaluate grad_f(x)
%     function [Hv] = Hessvec(v, x, usr_par)   evaluate Hess_f(x) * v
%
% and
%     function [usr_par] = xnew( x, iter, usr_par)
%     function [x1x2] = xprod( x1, x2, usr_par).
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
%             if newton_cg  is called with options =[], then defauls
%                   are used.
%         
%   usr_par problem specific information. usr_par is not referenced in
%           newton_cg, but passed to the user provided fucntions fval,
%           grad, ....
%
%
% output:
%    x      final approximate solution
%    iter   number of Newton iterations
%    iflag  error flag
%           iflag = 0, Successful  norm(g) < gtol
%           iflag = 1, linesearch returned with error
%           iflag = 2, CG did not compute a descent direction
%           iflag = 3, reached max number of iterations without achieving
%                      gradient norm > gtol
%           iflag = 4, step norm < stol, but gradient norm > gtol
%                 
%    
%   Matthias Heinkenschloss
%   Department of Computational and Applied Mathematics
%   Rice University
%   June 6, 2008
%
%
function [ x, iter, iflag ]  = newton_cg( x, options, usr_par)


%initialization
iflag       = 0;
iter        = 0;

%% set linear search parameters
iprint  = 0;
fid     = 1;
maxit   = 100;
gtol    = 1.e-8;
stol    = 1.e-8;
if( ~isempty(options) )
    if( isfield(options,'iprint') ); iprint = options.iprint; end
    if( isfield(options,'fid') );    fid    = options.fid; end
    if( isfield(options,'maxit') );  maxit  = options.maxit; end
    if( isfield(options,'gtol') );   gtol   = options.gtol; end
    if( isfield(options,'stol') );   stol   = options.stol; end
end


if( iprint > 0 )
   fprintf(1,' newton_cg \n')
   fprintf(1,['     k    f(x_k)    ||gradf(x_k)||    ||s_k||  ', ...
              '     stepsize     # CG iters  CG-flag\n'])
end

usr_par  = xnew( x, iter, usr_par);
f        = fval(x, usr_par);
g        = grad(x, usr_par);
gnorm    = sqrt(xprod( g, g, usr_par));
stepsize = 1;
snorm    = 2*stol;

%stopping criteria for newton
while ( ~(gnorm < gtol | stepsize*snorm < stol ) ...
        & iter <= maxit & iflag == 0)
   
   % Compute search direction
   tolcg    = min(gnorm^2, 0.01*gnorm);
   max_itcg = 2*length(x);
   s        = zeros(size(x));
   [s , res, itercg, iflagcg] = mycg( @(s)Hessvec(s,x,usr_par), s, -g, [], max_itcg, tolcg, usr_par); 
   
   % check CG error flag
   if( itercg == 0 & (iflagcg == 2 | iflagcg == 3) )
      % Negative eigenvalues was detected in Hessian. Take negative 
      % gradient direction if this occurred in first CG iteration
      s = -g;
   end
   % check descent property
   gs = xprod( g, s, usr_par );
   
   if( gs >= 0 )
       if( iprint > 0 )
          fprintf(fid,' Newton iteration %d \n', iter)
          fprintf(fid,' Inexact Newton direction s is not a descent direction')
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
   
   f_c      = f;  % save current function value for printing
   gnorm_c  = gnorm;
   stepsize = 1;  % initial stepsize
   [stepsize, x, f, iter_ls, iflag_ls, usr_par] = lnsrch_bt( x, s, gs, f, stepsize, ln_options, usr_par );  
%    [stepsize, x, f, iter_ls, iflag_ls, usr_par] = lnsrch_arm( x, s, gs, f, stepsize, ln_options, usr_par );  
   
   if( iprint > 1 )
   fprintf(fid,' newton_cg \n')
   fprintf(fid,['     k    f(x_k)    ||gradf(x_k)||    ||s_k||  ', ...
              '     stepsize     # CG iters  CG-flag\n'])
   end
   if( iprint > 0 )
      fprintf(fid,'  %4d  %12.6e  %12.6e  %12.6e  %12.6e    %6d    %6d \n', ...
                  iter, f_c, gnorm_c, snorm, stepsize, itercg, iflagcg )
   end
   
   if( iflag_ls > 0 )
      if( iprint > 0 )
          fprintf(fid,' Newton iteration %d \n', iter)
          fprintf(fid,' line search returned with flag = %d; terminate \n', iflag_ls)
      end
      iflag = 1;
      return
   end
   
   g        = grad(x, usr_par);
   gnorm    = sqrt(xprod( g, g, usr_par));
   

   iter  = iter+1;
end

if( iprint > 0 & iflag == 0 )
      fprintf(fid,'  %4d  %12.6e  %12.6e \n\n', iter, f, gnorm )
end

if ( gnorm > gtol & iter >= maxit)
   iflag = 3;
end
if ( gnorm > gtol & stepsize*snorm >= stol)
   iflag = 4;
end

   
   
   



