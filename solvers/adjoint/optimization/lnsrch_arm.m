function [stp, x, f, iter, iflag, usr_par] = lnsrch_arm( x, d, dg, f, stp, options, usr_par )
%
%   [stp, x, f, iter, iflag, usr_par] = lnsrch_bt( x, d, dg, f, stp, options, usr_par )
%
%   Use backtracking to find a step size  stp > 0  such that
%  
%          f(x+stp*d) <= f(x) + ftol * stp* < gradf(x), d >,
%  
%   where gradf(x) is the gradient at x and d is a search direction.
%  
%   This function calls the following usr provided functions
%     function [fx] = f(x, usr_par)        
%     function [usr_par] = xnew( x, iter, usr_par).
%
%   We use a imple ASrmijo line-search to find the step size.
%  
%   Inputs   
%      x            current iterate
%      d            search direction
%      dg           < gradf(x), d >
%      f            objective function value at x
%      stp          initial step size
%      options      line search parameters
%                      option.iprint   print if ~=0, default = 0
%                      option.fid      file identifier, default = 1
%                      option.stpmin   minimum step length, default = 1.e-8
%                      option.ftol     ftol in sufficient decrease, default = 1.e-4
%                   if lnsrch_bt  is called with options =[], then defauls
%                   are used.
%      usr_par      user parameters. These are never referenced
%                   in this function, but are passed to the
%                   user provided functions f and xnew.
%  
%   Outputs
%      stp          final step size
%      x            x + stp * d 
%      f            objective function at x+ stp * d
%      iter         number of iterations in lnsrch_arm
%                   (=number of fctn. evals. in lnsrch_arm) 
%      iflag        return flag
%                   iflag = 0 line search was successful
%                   iflag = 1 d is not a descent direction
%                   iflag = 2 initial stp < 0
%                   iflag = 3 maximum number of LS 
%                             iterations exceeded
%    
%      usr_par      user parameters. may be changed inside this function
%                   through a call to the user provided function xnew.
%
%
%   Matthias Heinkenschloss
%   Department of Computational and Applied Mathematics
%   Rice University
%   June 6, 2008
          

%% set linear search parameters
iprint  = 0;
fid     = 1;
stpmin  = 1.e-8;
ftol    = 1.e-4;
if( ~isempty(options) )
    if( isfield(options,'iprint') ); iprint = options.iprint; end
    if( isfield(options,'fid') );    fid    = options.fid; end
    if( isfield(options,'stpmin') ); stpmin = options.stpmin; end
    if( isfield(options,'ftol') );   ftol   = options.ftol; end
end

%%  check inputs
if( dg >= 0 )       %  d is not a descent direction
    iflag = 1;
    return
end
if( stp <= stpmin ) %  initial step too small
    iflag = 2;
    return
end

%%
if( iprint > 0 )
    fprintf(fid,' lnsrch_arm with \n')
    fprintf(fid,' option.stpmin = %12.5e \n', stpmin)
    fprintf(fid,' option.ftol   = %12.5e \n', ftol)
    fprintf(fid,' iter    stepsize       Obj_fctn       required decrease \n'); 
end

%% initialize
iter   = 0;
f0     = f;
x0     = x;

%% start linear search loop
while (stp > stpmin) 

    %  trial iterate 
    x = x0 + stp * d;

    %  compute the function value at x+stp*d
    usr_par  = xnew( x, iter, usr_par);
    f        = fval( x, usr_par);
   
    %  print iteration info
    if( iprint > 0 )
        fprintf(fid, '%4d    %12.5e    %12.5e   %12.5e \n', ...
                      iter, stp, f, f0 + stp*ftol*dg )
    end

    %  check decrease condition and return if successful
    if( f < f0 + stp*ftol*dg )
        iflag  = 0;
        return
    end

    % reduce step size 
    stp  = 0.5*stp;
    iter = iter + 1;
end

%% stp < stpmin
iflag  = 3;


