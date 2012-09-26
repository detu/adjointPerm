function z = blackoilUpwFE(z, tf, q, F, pv, fluid, PVTTAB, varargin)
%Explicit single point upwind transport solver for Black Oil flow.
%
% SYNOPSIS:
%    z = blackoilUpwFE(z, t, q, gm, pv, fluid, table)
%    z = blackoilUpwFE(z, t, q, gm, pv, fluid, table, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function blackoilUpwFE computes the solution to the Black-Oil transport
%   equation
%
%        z_t + f(s)_x/B = q, with z = s / B,
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,t] while honouring a CFL time step restriction.
%
% PARAMETERS:
%   z       - An N-by-3 array (matrix) of initial phase densities (at
%             time=0).  Here, N is the number of cells in the model and the
%             entries in 'z' are interpreted as
%
%                z(c,1) <-> Density of aqueous phase (water) in cell c
%                z(c,2) <-> Density of liquid phase  (oil)   in cell c
%                z(c,3) <-> Density of gaseous phase (gas)   in cell c
%
%   t       - End point of time integration interval (i.e., final time)
%
%   q       - Accumulated sources (e.g., contributions from wells and
%             boundary conditions).  One scalar value for each cell in the
%             model.
%
%   gm      - Three-element cell array of upwind/inflow matrices of
%             phase-dependent fluxes into each cell, created (e.g.) by
%             function INFLOW_BO.  One matrix for each of the phases 'Aqua'
%             (1), 'Liquid' (2), and 'Vapor' (3).
%
%             Specifically, the entry gm{p}(i,j) is the (positive) phase
%             flux of phase 'p' from cell j to cell i.
%
%   pv      - Reservoir pore volumes.  One scalar value for each cell in the
%             model.
%
%   fluid   - Data structure describing the fluids in the problem.  Typically
%             created by function EVALPVT.
%
%   table   - PVT table containing (raw) tabulated PVT data.
%             Must be suitable for passing to function getRelPerm.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - dt      -- Time step size.  Double.  Default value = t.
%               - verbose -- Whether or not to emit progress reports during
%                            process.  Logical.  Default value: False.
%
% RETURNS:
%   z       - Phase densities for each cell in the model at time=t.
%
% EXAMPLE:
%   blackOilTube
%
% SEE ALSO:
%   inflow_bo, getRelPerm, blackOilUpwBE.

%{
#COPYRIGHT#
%}

% $Id: blackoilUpwFE.m 1953 2009-03-31 10:54:12Z bska $

opt = struct('dt', tf, 'verbose', false);
opt = merge_options(opt, varargin{:});

dt = opt.dt;
nc = size(F{1}, 1);
np = numel(F);      % Number of phases

for p = 1 : np,
   d    = sum(F{p}).';   % d(c) == amount of phase 'p' leaving cell 'c'.
   F{p} = F{p} - spdiags(d, 0, nc, nc);
end


if opt.verbose,
   fprintf('Solving transport using %d time steps...', ceil(tf / dt));
end


f  = zeros([nc, np]);
s  = zeros([nc, np]);
dz = zeros([nc, np]);
H  = @(f,dz,dt) bsxfun(@times, dt ./ pv, ...
                       dz - max(q,0) - bsxfun(@times, min(q,0), f));

t = 0;
wheel  = ['-', '\', '|', '/'];
iwheel = 0;
while t < tf,
   iwheel = mod(iwheel, 4) + 1;    % == 1, 2, 3, 4, 1, 2, 3, 4, ...
   fprintf('\b%s', wheel(iwheel));

   % Compute saturations for purpose of relative permeabilities.
   %
   s(:) = fluid.B .* z;
   s(:) = bsxfun(@rdivide, s, sum(s,2));

   % Compute fractional flow functions to determine the flow of the
   % individual phases.
   %
   f(:) = getRelPerm(s, PVTTAB);
   f(:) = bsxfun(@rdivide, f, fluid.mu);
   f(:) = bsxfun(@rdivide, f, sum(f,2));

   % Actual transport of masses during (short) timestep dt according to
   % phase fluxes.
   %
   for p = 1 : np, dz(:,p) = F{p} * f(:,p); end

   z(:) = z - H(f, dz, dt);

   t  = t + dt;
   dt = min(dt, tf - t);
end
fprintf('\b.\n');
