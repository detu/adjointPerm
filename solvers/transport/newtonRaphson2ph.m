function resSol = newtonRaphson2ph(resSol, tf, F, Jac, update, opt)
%Solve non-linear equation F(s)=0 using Newton-Raphson method.
%
% SYNOPSIS:
%   s = newtonRaphson2ph(resSol, tf, F, Jac, update, opt)
%
% DESCRIPTION:
%   Solve general nonlinear equation F(s)=0 using a Newton-Raphson
%   iteration with user-specified UPDATE scheme to adjust iterate.
%   Specifically, the k-th iteration is defined schematically as
%
%      ds^k = J \ F
%      s^k  = update(s^k-1, ds^k, ...)
%
% PARAMETERS:
%   resSol      -
%   tf     -
%   F      -
%   Jac    -
%   update -
%   opt    -
%
% RETURNS:
%   resSol -

%{
#COPYRIGHT#
%}

% $Date: 2009-10-09 11:37:04 +0200 (fr, 09 okt 2009) $
% $Revision: 2969 $

   report = struct('success', true,...
                   'iterations', 0, ...
                   'sub_steps', 0, ...
                   'failed_steps', 0);
                   

   % maxnewt?, nltol?
   mints   = pow2(tf, -opt.tsref);
   dispif(opt.verbose, '\n\n');
   dispif(opt.verbose, [repmat('-', [1, 70]), '\n']);
   dispif(opt.verbose, '  Time interval            iter          residual \n');
   dispif(opt.verbose, [repmat('-', [1, 70]), '\n']);

   [t, dt, dtprev, count] = deal(0.0, tf, -tf, 0);
   
   %-----------------------------------------------------------------------
   % Main time loop (ideally, just a single step) -------------------------
   %
   while t < tf && dt > mints,
      dt = min(dt, tf - t);
      %setDirections(resSol);
      
      %--------------------------------------------------------------------
      % Outer controlling loop (sub step size &c) -------------------------
      %
      redo_newton = true;
      while redo_newton,
         dispif(opt.verbose, '[%f, %f]:', t / tf, (t + dt) / tf); nc = 0;
         s0 = resSol.s(:,1);
         sn = resSol; sn.s(:,1) = 0.5;

         %-----------------------------------------------------------------
         % Initialise NR algorithm ----------------------------------------
         %
         res     = F(sn, s0, dt);
         err     = norm(res, inf);
         nwtfail = err > opt.nltol;
         linfail = false;
         it      = 0;

         %-----------------------------------------------------------------
         % Execute inner NR algorithm -------------------------------------
         %
         while nwtfail && ~linfail && it < opt.maxnewt,
            ds  = - opt.LinSolve(Jac(sn, dt), res);
            [sn, res, alph, linfail] = update(sn, s0, ds, dt, err);
            it      = it + 1;
            err     = norm(res, inf);
            nwtfail = err > opt.nltol;

            dispif(opt.verbose, repmat('\b', [1, nc]));
            nc = dispif(opt.verbose, '\t%4d\t %2.2f\t %5.5e', it, alph, err);
            %dispif(opt.verbose, '\n\t%4d\t %2.2f\t %5.5e', it, alph, err);          
            report.iterations = report.iterations +1;
         end

         %-----------------------------------------------------------------
         % Determine if NR succeeded or not -------------------------------
         %
         count = count - 1;
         if nwtfail,
            dispif(opt.verbose, '\tReducing step.');

            if count > 0 && dtprev > 0,
               dt = dtprev;
            else
               dt = dt / (1.5 + 0.5);
            end
            count = 5;
            report.failed_steps = report.failed_steps + 1;
         else
            redo_newton = false;
            t = t + dt;  dtprev = -dt; %invalid dtprev.
            if it == 0,
               dispif(opt.verbose, '\t%4d\t %s\t %5.5e', it, '-', err);
               dispif(opt.verbose, '\t NB: err <= ntol.');
            end
            if t<tf && it <= 5  && count < 0, % arbitrary threshold
               dispif(opt.verbose, '\tIncreasing step.');

               dtprev = dt;
               dt     = min((1 + 0.5) * dt, tf - t);
               count  = 5;
            end
         end
         dispif(opt.verbose, '\n');
      end

      % Time step [t, t+dt] was successful
      resSol = sn;
      report.sub_steps = report.sub_steps+1;
   end
   
   report.success = ~(t<tf);
   if ~report.success,
      resSol.s(:,1) = NaN(size(resSol.s, 1), 1);
   end
   
   % Update oil saturation if applicable.
   if size(resSol.s,2) > 1,
      resSol.s(:,2) = 1 - resSol.s(:,1);
   end

   if opt.verbose, dispReport(report, dt, mints); end
end




function dispReport(r, dt, mints)
   fprintf('\n');
   fprintf([repmat('-', [1, 70]), '\n']);
   fprintf(['Transport step used %d iterations in %d sub steps ',...
            '(%d failed steps)\n'], r.iterations, r.sub_steps,   ...
           r.failed_steps);
   if ~r.success,
      fprintf(['Transport step FAILED due to timestep %g < %g.\n'], ...
              dt, mints);
             
   end   
end
