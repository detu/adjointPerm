function z = newtonRaphson(z, tf, F, Jac, update, opt)
%Solve non-linear equation F(z)=0 using Newton-Raphson method.
%
% SYNOPSIS:
%   z = newtonRaphson(z, tf, F, Jac, update, opt)
%
% DESCRIPTION:
%   Solve general nonlinear equation F(z)=0 using a Newton-Raphson
%   iteration with user-specified UPDATE scheme to adjust iterate.
%   Specifically, the k-th iteration is defined schematically as
%
%      dz^k = J \ F
%      z^k  = update(z^k-1, dz^k, ...)
%
% PARAMETERS:
%   z      -
%   tf     -
%   F      -
%   Jac    -
%   update -
%   opt    -
%
% RETURNS:
%   z -

%{
#COPYRIGHT#
%}

% $Id: newtonRaphson.m 1953 2009-03-31 10:54:12Z bska $

   mints   = pow2(tf, -opt.tsref);
   dispif(opt.verbose, '\n\n');
   dispif(opt.verbose, [repmat('-', [1, 70]), '\n']);
   dispif(opt.verbose, '  Time interval                 iter    residual \n');
   dispif(opt.verbose, [repmat('-', [1, 70]), '\n']);

   [t, dt, count] = deal(0.0, tf, 0);

   %-----------------------------------------------------------------------
   % Main time loop (ideally, just a single step) -------------------------
   %
   tic
   while t < tf && dt > mints,
      dt = min(dt, tf - t);

      %--------------------------------------------------------------------
      % Outer controlling loop (sub step size &c) -------------------------
      %
      redo_newton = true;
      while redo_newton,
         dispif(opt.verbose, '[%12.4f, %12.4f]:', t, t + dt); nc = 0;

         [z0, zn] = deal(z);

         %-----------------------------------------------------------------
         % Initialise NR algorithm ----------------------------------------
         %
         res     = F(zn, z0, dt);
         err     = norm(res, inf);
         nwtfail = err > opt.nltol;
         linfail = false;
         it      = 0;

         %-----------------------------------------------------------------
         % Execute inner NR algorithm -------------------------------------
         %
         while nwtfail && ~linfail && it < opt.maxnewt,
             J   = Jac(zn, dt);
             dz  = - J \ res;
             [zn, res, linfail] = update(zn, z0, dz, dt, err);
             %zn  = max(0, zn + reshape(dz, size(zn'))');
             zn  = max(zn, 0);
             res     = F(zn, z0, dt);
             it      = it + 1;
             err     = norm(res, inf);
             nwtfail = err > opt.nltol;

             dispif(opt.verbose, repmat('\b', [1, nc]));
             nc = dispif(opt.verbose, '  %4d    %5.5e', it, err);
         end

         %-----------------------------------------------------------------
         % Determine if NR succeeded or not -------------------------------
         %
         count = count - 1;
         if nwtfail,
            dispif(opt.verbose, '\tReducing step.');

            dt    = dt / 2;
            count = 5;
         else
            redo_newton = false;
            t = t + dt;
            if t<tf && it <= 5  && count < 0, % arbitrary threshold
               dispif(opt.verbose, '\tIncreasing step.');

               dt = min(1.5 * dt, tf - t);
            end
         end
         dispif(opt.verbose, '\n');
      end

      % Time step [t, t+dt] was successful
      z = zn;
   end
   tocif(opt.verbose)

   if (t < tf),
      dispif(opt.verbose, 'Unable to integrate to final time, t=%f\n', tf);
      z(:) = NaN;
   end
end
