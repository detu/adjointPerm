function fluid = initCompressibleFluid(rho0, p0, c0, mu0, varargin)
%Construct PVT and relperm functions, constant compr., single phase.
%
% SYNOPSIS:
%   fluid = initCompressibleFluid(rho0, p0, c0, mu0)
%
% PARAMETERS:
%   rho0 - Fluid reference density, in units of kg/m^3.
%
%   p0   - Pressure, in units of Pascal, at which fluid reference
%          parameters (density, compressibility, and viscosity) are
%          defined.
%
%   c0   - Fluid reference compressibility at pressure 'p0'.
%
%   mu0  - Fluid reference viscosity, in units of Pa s, at pressure 'p0'.
%
% RETURNS:
%   fluid - A 1-by-1 structure array of PVT and relative permeability
%           evaluation functions.  The individual evaluation functions
%           (structure fields) are,
%             - pvt -- Evaluates pvt data for the single phase compressible
%                      fluid.  Specifically, given an n-by-1 vector of
%                      pressures p and an n-by-1 vector of surface volumes
%                      z, the call
%
%                         [c, rho, mu, u] = fluid.pvt(p, z)
%
%                      computes, respectively, n-by-1 phase
%                      compressibilities (c), densities (rho), viscosities
%                      (mu), and reservoir volumes (u).
%
%              - relperm --
%                      Evaluates one-phase relative permeability curves.
%                      Specifically, given an n-by-1 vector of fluid
%                      saturations, s (implicitly assumed s = ones([n,1])),
%                      the call
%
%                         kr = fluid.relperm(s)
%
%                      computes an n-by-1 array of relative permeability
%                      values.  As 'fluid' is a single phase fluid model,
%                      all relative permeability values are one, i.e.,
%
%                         kr = ones(size(s))
%
% EXAMPLE:
%
% SEE ALSO:
%   initBlackoilFluid.

%{
#COPYRIGHT#
%}

% $Id: initCompressibleFluid.m 2926 2009-10-01 14:38:31Z bska $

   function [c, rho, mu, u, R] = pvtfun(varargin)
      % Poor man's class replacement
      if nargin == 1,
         str = varargin{1};
         assert(ischar(str));

         if strcmpi(str, 'info'),
            assert (nargout == 0);
            disp('single compressible fluid');
            return;
         elseif strcmpi(str, 'surfacedensity'),
            assert (nargout <= 1);
            c = rho0;
            return
         elseif strcmpi(str, 'miscible'),
            assert (nargout <= 1);
            c = false;
            return;
         elseif strcmpi(str, 'names'),
            c = '';
            return
         else
            error('unsupported mode');
         end
      elseif nargin == 2 && isnumeric([varargin{1:2}]),
         p = varargin{1};
         z = varargin{2};
      else
         error('Unsupported input.');
      end

      if numel(p) ~= size(z, 1),
         error('There must be one pressure for each row in mass table');
      end

      c    = repmat(c0,  [numel(p), 1]);
      rho  = rho0 .* exp(c .* (p - p0));
      mu   = repmat(mu0, [numel(p), 1]);
      u    = rho0 .* z ./ rho;
      R    = zeros([numel(p), 1]);
   end

   fluid.pvt     = @pvtfun;
   fluid.relperm = @(s) ones([size(s,1), 1]);
end
