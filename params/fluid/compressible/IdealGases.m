classdef IdealGases < CompressibleFluidBase
   % Simple model of an ideal gas.
   % The ideal gas is assumed to follow the rule: pV = 1
   % Input parameters: viscosity, surface density, and surface pressure.
   % 

   % $Date: 2009-10-05 11:10:25 +0200 (ma, 05 okt 2009) $
   % $Revision: 2945 $

   properties (GetAccess=public, SetAccess=protected)
      surfaceDensity  = [];
      surfacePressure = [];
      info            = [];
      viscosity       = [];
      miscible        = [];
   end

   methods (Access=public)
      function self = IdealGases(mu, rhoS, pS)
         self.viscosity       = mu;
         self.surfaceDensity  = rhoS;
         self.surfacePressure = pS;
         self.info            = ...
            sprintf(['Ideal gas(es):\n'                        ...
                     'Reference pressure: %g (Pa)\n',      ...
                     'At reference pressure:\n',           ...
                     '  Density:         [%6g] (Kg/m3)\n', ...
                     '  viscosity:       [%6g] (Pa·s)\n'], ...
                     self.surfacePressure, self.surfaceDensity(:), ...
                     self.viscosity(:));

         self.miscible = false([1, numel(rhoS)]);
      end

      function [kr, dkr] = relperm(self, s)
         assert (size(s,2) == numel(self.surfaceDensity),  ...
                 ['Saturation must have the same number of ', ...
                  'components as the density'])
         kr  = s;
         dkr = ones(size(s));
      end

      function [c, rho, mu, u, frac] = pvt(self, p, z)
         n    = numel(p);
         rho  = bsxfun(@times, self.surfaceDensity, p/self.surfacePressure);
         c    = bsxfun(@rdivide, ones([1, numel(self.surfaceDensity)]), p);
         mu   = repmat(self.viscosity,       [n,1]);
         u    = bsxfun(@rdivide, bsxfun(@times, z, self.surfaceDensity), rho);
         frac = zeros(size(z));
      end
   end
end
