function resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, ...
                                    pv, fluid,  varargin)
%Explicit single point upwind solver for two-phase flow, including gravity.
%
% SYNOPSIS:
%   resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, porvol,
%                              fluid)
%   resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, porvol,
%                              fluid, G, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwFEGrav solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf].
%
%   The upwind forward Euler discretisation of the Buckley-Leverett model
%   can be written as:
%
%     s^(n+1) = s^n - (dt./pv)*((H(s^n) - max(q,0) - min(q,0)*f(s^n))
%
%   where
%        H(s) = (flux + grav*diag(A_o*lam_o(s))
%               *A_w*lam_w(s)./(A_w*lam_w(s)+A_o*lam_o(s)),
%
%   pv is the porevolume, lam_l is the mobility for face l, f is
%   Buckely-Leverett fractional flow function, while A_o and A_w are index
%   matrices that determine the upstream mobility.
%
% PARAMETERS:
%   resSol  - Reservoir solution structure containing valid water
%             saturation resSol.s(:,1) with one value for each cell
%             in the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   q       - Accumulated sources (e.g., contributions from wells and
%             boundary conditions).  One scalar value for each cell in the
%             model.  The source rates are assumed to be measured in units
%             of m^3/s.
%
%   flux    - Inflow matrix of fluxes into each cell, created
%             by function initTransport.  Size
%             G.cells.num-by-(G.faces.num-boundary faces)
%
%             Specifically, the entry flux(i,j) is the flux over
%             (active) face j from cell i.  The fluxes are assumed to be
%             measured in units of m^3/s.
%
%   grav    - Matrix with gravity contribution for each face, created by
%             function initTransport.  Same size as flux.
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
%                           seconds, for instance determined by function
%                           initTransport.  Default value = tf.
%                           NB: The explicit scheme is only stable provided
%                           that dt satisfies a CFL time step restriction.
%
%               - satWarn - Tolerance level for saturation warning.
%                            Default value: satWarn = sqrt(eps).
% RETURNS:
%   resSol - Reservoir solution with updated saturations, resSol.s.
%
% SEE ALSO:
%   twophaseUpwFE, initTransport, explicitTransport, implicitTransport.

% TODO:
%   - implement gravity effects for pressure boundary and wells

%{
#COPYRIGHT#
%}

% $Id: twophaseUpwFEGrav.m 2963 2009-10-09 05:06:18Z jrn $

assert (size(resSol.s,2)<3 || all(resSol.s(:,3)==0));

% Verify that the flux matrix is correctly sized.  The output from
% function initTransport depends on presence or absence of gravity effects.
%
if size(flux) ~= size(grav),
   error('twophaseUpwFEGrav:InputSize:flux',                    ...
        ['Use twophaseUpwFEGrav only in presence of gravity, ', ...
         'else use twophaseUpwFE.'])
end

opt = struct('dt', tf, 'verbose', false, 'satWarn', sqrt(eps));
opt = merge_options(opt, varargin{:});
dt  = opt.dt;
dt  = min(dt, tf);

assert (dt > 0 && tf > 0);

[v_darcy, neighbors, normals, upw_inx, g_vec] = initFaceMob(G, resSol, ...
                                                           flux, grav);

H = @(f,dz,dt) ((dt ./ pv).*(dz - max(q,0) - min(q,0).*f));
t = 0;

if opt.verbose,
%   h = waitbar(0, ['Solving transport using ' ...
%                   num2str(ceil(tf/dt)), ' time steps...']);
end
while t < tf,
   % Compute cell mobility and fractional flow.
   mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
   f_w = mob(:,1)./sum(mob,2);

   % Initialize face mobility.
   if nnz(grav)
      [iw, io] = findFaceMobIx(upw_inx, mob, g_vec, ...
                                  v_darcy, neighbors, normals, fluid);
      faceMob = [mob(iw,1) mob(io,2)];
   else
      faceMob = mob(upw_inx, :);
   end

   fw_face = faceMob(:,1) ./ sum(faceMob,2);
   %Remove possible NaNs
   fw_face(sum(faceMob,2)==0) = 0;

   dz = flux*fw_face + grav*(faceMob(:,2) .* fw_face);

   % Explicit computation of saturation
   resSol.s(:,1) = resSol.s(:,1) - H(f_w, dz, dt);

   % Correct possible bad saturations
   if opt.verbose,
      %waitbar(t/tf,h);
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

%if opt.verbose, close(h), end

if size(resSol.s,2) > 1,
   % Update oil saturation:
   resSol.s(:,2) = 1 - resSol.s(:,1);
end
