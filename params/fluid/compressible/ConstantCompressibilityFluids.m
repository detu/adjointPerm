classdef ConstantCompressibilityFluids < CompressibleFluidBase
   %Class representing compressible fluids with constant compressibility.

   % $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
   % $Revision: 2926 $

   properties(GetAccess=public, SetAccess=protected)
      p_ref           = [];
      compressibility = []

      surfaceDensity  = [];
      info            = [];
      viscosity       = [];
      miscible        = [];
   end


   methods(Access=public)
      function self = ConstantCompressibilityFluids(c, rho, p, varargin)
         assert (numel(p) == 1);
         self.p_ref           = p;
         self.compressibility = reshape(c, 1, []);
         self.surfaceDensity  = reshape(rho, 1, []);
         if nargin > 3,
            self.viscosity = reshape(varargin{1}, 1, []);
         else
            self.viscosity = ones([1, numel(c)]) * centi * poise();
         end
         self.info = self.getInfo();
         self.miscible = false([1, numel(c)]);
      end

      function kr = relperm(self, s)%#ok
         %assert(size(s,2)==1,  'Saturation must be a column vector')
         kr  = s.^2;
      end

      function [c, rho, mu, u, frac] = pvt(self, p, z)
         n    = numel(p);
         t    = bsxfun(@times, p(:)-self.p_ref, self.compressibility);
         rho  = bsxfun(@times, self.surfaceDensity, exp(t));
         c    = repmat(self.compressibility, [n,1]);
         mu   = repmat(self.viscosity,       [n,1]);
         u    = bsxfun(@rdivide, bsxfun(@times, z, self.surfaceDensity), rho);
         frac = zeros(size(z));
      end
   end


   methods (Access=protected)
      function s = getInfo(self)
         n = numel(self.surfaceDensity);
         fmt = ['  %-15s: [', repmat('%6g ', [1, n]), '\b] %s\n'];
         s = sprintf(['Constant compressibility fluid(s):\n',...
                      'Reference pressure: %g (Pa)\n'], self.p_ref);
         s = [s, sprintf('At reference pressure:\n'), ...
              sprintf(fmt, 'density', self.surfaceDensity(:), 'kg/m3'), ...
              sprintf(fmt, 'compressibility', self.compressibility(:), '1/Pa'), ...
              sprintf(fmt, 'viscosity',  self.viscosity(:), 'Pa·s')];
      end
   end
end
