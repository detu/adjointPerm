
function [x, res, iter, flag] = mycg( A, x, b, M, max_it, tol, usr_par ) 
%
%  CG:
%     computes an approximate solution of   min 0.5*<Ax,x> - <b,x> 
%     using the preconditioned conjugate gradient algorithm.
%
%     The user has provide a functions 
%               function [x1x2] = xprod( x1, x2, usr_par)
%     which evaluates the inner product <x1,x2>  of x1 and x2.
%     
%
%  Usage:
%
%    [x, res, iter, flag] = mycg( A, x, b, M, max_it, tol )
%
%
%  Input parameters
%
%  A       the system matrix A.
%          is either a real(n,n) symmetric positive definite matrix or 
%          the name of a function that computes y = A * x.
%          The calling sequence of this function is [y] = A(x,usr_par).
%          The operator A has t be symmetric, i.e, it has to satisfy
%          <A*x1, x2> = <x1, A*x2> for all x1, x2
%
%  x       real(n)
%          the initial approximate solution.
%
%  b       real(n)
%          the right-hand side. 
%
%  M       the preconditioner
%          is either a symmetric positive definite n by n matrix or the 
%          name of a function that solves M y = x.
%          The calling sequence of this function is  [y] = M(x,usr_par).
%          No preconditioning is applied when M = [].
%          The operator M has to be symmetric, i.e, it has to satisfy
%          <M*x1, x2> = <x1, M*x2> for all x1, x2
%
%  max_it  integer
%          maximum number of iterations.
%
%  tol     real
%          the stopping tolerance. 
%          the algorithm attempts to generate x such that   
%          r^T M^(-1) r <= tol,  where r = b-Ax.
%
%   usr_par problem specific information. usr_par is not referencced in
%           mycg, but passed to the function evaluating A, M,...
%
%  Output parameters: 
%
%  x       Final approximate solution. 
%
%  flag    integer
%          Error flag, as follows: 
%              0 => normal termination; desired tolerance met. 
%              1 => desired tolerance not met after max_it iterations.
%              2 => in iteration iter an inner product of the form
%                   v' M^(-1) v was not positive, so the preconditioning 
%                   matrix M  does not appear to be positive definite.
%                   x  will not contain an acceptable solution.
%              3 => in iteration iter an inner product of the form
%                   <q, A q> was not positive, so the matrix A  does not 
%                   appear to be positive definite.
%                   x  will not contain an acceptable solution.
%
%  res     real
%          vector of residual norms. 
%          res(i+1) is the preconditioned residual norm of iteration i.
%          res(i+1) = sqrt( <r, M^(-1) r> )
%
%  iter    integer
%          contains the total number of iterations in CG.
%          the number of  matrix vector products is given by  iter + 1.
% 
%
%   Matthias Heinkenschloss
%   Department of Computational and Applied Mathematics
%   Rice University
%   June 6, 2008
%   
%

  % initialization
  flag = 0;        

  % compute initial residual
  if( isstr(A) )
        r = b - feval( A, x, usr_par );
  elseif( isa(A, 'function_handle') )
        r = b - A(x);
  else
        r = b - A*x;
  end
  %res(iter+1) = norm( r );     

  % begin iteration
  for iter = 0:max_it 

     % apply  preconditioner
     if( isstr(M) )
        z = feval( M, r, usr_par );
     elseif( isa(M, 'function_handle') )
        z = M(r);
     elseif( ~isempty(M) )
        z = M\r;
     else
        z = r;
     end
     rho = xprod( r, z, usr_par );
     if( rho < 0 ); flag = 2; return; end
     
     res(iter+1) = sqrt(rho); 
     
     % check convergence  
     if ( res(iter+1) <= tol ); return; end

     % direction vector
     if ( iter > 0 ),      
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

     if( isstr(A) )
         q = feval( A, p, usr_par );
     elseif( isa(A, 'function_handle') )
         q = A(p);
     else
         q = A*p;
     end
     alpha = rho / xprod( p, q, usr_par );
     
     %disp([rho,xprod( p, q, usr_par ),alpha])
     if( alpha < 0 ); 
         % negative curvature detected. Return iterate computed thus far
         flag = 3; 
         return;
     end
     
     % update approximation vector
     x = x + alpha * p;             

     % compute residual
     r = r - alpha*q;       

     rho_1 = rho;

  end
  
  % no convergence
  flag = 1;       

% END mycg.m
