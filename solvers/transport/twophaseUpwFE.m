function resSol = twophaseUpwFE(resSol, tf, q, F, pv, fluid, varargin)
%Explicit single point upwind solver for two-phase flow, no gravity.
%
% SYNOPSIS:
%   s = twophaseUpwFE(resSol, tf, q, F, porvol, fluid)
%   s = twophaseUpwFE(resSol, tf, q, F, porvol, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwFE solves the Buckley-Leverett transport equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf]. The stability of the scheme depends on whether the
%   internal timestep satiesfies a CFL time step restriction.
%
% PARAMETERS:
%   resSol  - Reservoir solution structure containing valid (water)
%             saturation resSol.s(:,1) with one value for each cell in the
%             grid.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   q       - Accumulated sources (e.g., contributions from wells and
%             boundary conditions).  One scalar value for each cell in the
%             model.  The source rates are assumed to be measured in units
%             of m^3/s.
%
%   F       - Upwind/inflow matrix of fluxes into each cell, created (e.g.)
%             by function initTransport.  Specifically, the entry F(i,j) is
%             the (positive) flux from cell j to cell i.  The fluxes are
%             assumed to be measured in units of m^3/s.
%
%   porvol  - Reservoir pore volumes, measured in units of m^3.  One scalar
%             value for each cell in the model.
%
%   fluid   - Data structure describing the fluids in the problem.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%               - verbose - Whether or not time integration progress should
%                           be reported to the screen.
%                           Default value: verbose = false.
%
%               - dt      - Internal timesteps, measured in units of
%                           seconds. Default value = tf.
%                           NB: The explicit scheme is only stable provided
%                           that dt satiesfies a CFL time step restriction.
%
%               - satWarn - Tolerance level for saturation warning.
%                           Default value: satWarn = sqrt(eps).
%
% RETURNS:
%   resSol  - Reservoir solution structure with updated resSol.s(:,1).
%
% SEE ALSO:
%   twophaseUpwFEGrav, initTransport, explicitTransport, implicitTransport.

%{
#COPYRIGHT#
%}

% $Id: twophaseUpwFE.m 2044 2009-04-20 17:18:53Z bska $

% Verify that the flux matrix is correctly sized.  The output from
% function initTransport depends on presence or absence of gravity effects.
%
if size(F,1) ~= size(F,2),
   error('twophaseUpwFE:InputSize:F',                      ...
        ['Use twophaseUpwFE only in absence of gravity, ', ...
         'else use twophaseUpwFEGrav.'])
end
opt = struct('dt', tf, 'verbose', false, 'satWarn', sqrt(eps));
opt = merge_options(opt, varargin{:});

dt = opt.dt;

assert (dt > 0 && tf > 0);

nc = size(F, 1);
d  = sum(F).';   % d(c) == amount of fluid leaving cell 'c'.
F  = F - spdiags(d, 0, nc, nc);

H = @(f,dz,dt) ((dt ./ pv) .* (dz - max(q,0) - min(q,0).*f));

t = 0;

if opt.verbose,
   h = waitbar(0, ['Solving transport using ' ...
                   num2str(ceil(tf/dt)), ' time steps...']);
end
while t < tf,
   mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
   f_w = mob(:,1) ./ sum(mob,2);

   % Actual transport of masses during (short) timestep dt according to
   % phase fluxes.
   %
   dz = F * f_w;

   % Explicit computation of saturation.
   resSol.s(:,1) = resSol.s(:,1) - H(f_w, dz, dt);

   % Correct possible bad saturations
   if opt.verbose,
      waitbar(t/tf,h);
      sMax = max(resSol.s(:,1));
      if sMax> 1 + opt.satWarn,
         inx = find(resSol.s(:,1) > 1+opt.satWarn);
         disp('Saturation exceeds 1 in cells from:')
         fprintf(1,'%d %g\n', [inx, resSol.s(inx,1)]');
      end
      sMin = min(resSol.s(:,1));
      if sMin < -opt.satWarn,
         inx = find(resSol.s(:) < -opt.satWarn);
         disp('Saturation less than 0 in cells from:')
         fprintf(1,'%d %g\n', [inx, resSol.s(inx,1)]');
      end
   end
   resSol.s(resSol.s(:,1)>1,1) = 1 - eps;
   resSol.s(resSol.s(:,1)<0,1) = 0;

   t  = t + dt;
   dt = min(dt, tf - t);
end

if opt.verbose, close(h), end

if size(resSol.s,2) > 1,
   % Update oil saturation:
   resSol.s(:,2) = 1 - resSol.s(:,1);
end
