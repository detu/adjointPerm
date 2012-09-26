classdef Fluid < CompressibleFluidBase
   %Class representing a fluid composition as a collection of fluids.

   % $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
   % $Revision: 2926 $

   properties (GetAccess=public, SetAccess=protected)
      surfaceDensity = [];
      info           = [];
      viscosity      = [];
      miscible       = [];
   end

   properties (Access=protected)
      fluids         = {};
   end

   methods (Access=public)
      function self = Fluid(varargin)
         %Construct composite fluid model from individual components.
         %
         % SYNOPSIS:
         %   fluid = Fluid(f1, f2, ..., fn)
         %
         % PARAMETERS:
         %   f1, ..., fn - Individual fluids which collectively constitute
         %                 the complete fluid model in a simulation.
         %
         % RETURNS:
         %   fluid - Fully assembled fluid model.
         %
         % CAVEAT:
         %   At least one fluid component must be specified.

         assert (numel(varargin) > 0);
         self.fluids = varargin;
      end

      function varargout = relperm(self, s)
         %Compute relative permeability functions.
         %
         % SYNOPSIS:
         %    kr       = fluid.relperm(s)
         %   [kr, dkr] = fluid.relperm(s)
         %
         % PARAMETERS:
         %   s   - n-by-np array of fluid saturations.  One column for each
         %         of the 'np' fluid phases in the model.  Specifically,
         %         s(i,j) is the saturation of the 'j'th fluid component in
         %         the, typically, 'i'th grid cell.
         %
         % RETURNS:
         %   kr  - An n-by-np array of evaluated relative permeability
         %         functions.  Specifically, kr(i,j) is the relative
         %         permeability of fluid component 'j' evaluated at
         %         saturation s(i,j).
         %
         %   dkr - An n-by-np^2 array of numerically evaluated relperm
         %         derivatives.  Specifically, dkr(i, (j-1)*np+k) is the
         %         partial derivative of the 'j'th relative permeability
         %         curve with respect to the 'k'th fluid component all
         %         evaluated at saturation s(i,j).
         %         OPTIONAL.  Only returned if specifically requested.

         sz = numel(self.surfaceDensity);
         if size(s, 2) ~= sz,
            error('Expected %d saturation components, got %d.', ...
                  sz, size(s,2));
         end

         i = 0;
         a = cell([1, max(1, nargout)]);
         for f = 1 : numel(self.fluids),
            sz = numel(self.fluids{f}.surfaceDensity);
            ix = i + 1 : i + sz;

            [a{:}] = self.fluids{f}.relperm(s(:, ix));

            for k = 1 : numel(a), varargout{k}(:,ix) = a{k}; end

            i = i + sz;
         end
      end

      function varargout = pvt(self, p, z)
         %Compute pressure and mass-dependent derived fluid parameters.
         %
         % SYNOPSIS:
         %    c                    = fluid.pvt(p, z)
         %   [c, rho]              = fluid.pvt(p, z)
         %   [c, rho, mu]          = fluid.pvt(p, z)
         %   [c, rho, mu, u]       = fluid.pvt(p, z)
         %   [c, rho, mu, u, frac] = fluid.pvt(p, z)
         %
         % PARAMETERS:
         %   p -
         %   z - n-by-np array of fluid surface volumes.  One column for
         %       each of the 'np' fluid phases in the model.  Specifically,
         %       z(i,j) is the surface volume of the 'j'th fluid component
         %       in the, typically, 'i'th grid cell.
         %
         % RETURNS:
         %   c    - Array, n-by-np, of fluid phase compressibilities.
         %   rho  - Array, n-by-np, of fluid phase densities.
         %   mu   - Array, n-by-np, of fluid viscosities.
         %   u    - Array, n-by-np, of fluid reservoir volumes.
         %   frac - Review Battlestar Galactica.
         %
         % SEE ALSO:
         %   solveBlackOilWellSystem.

         sz = numel(self.surfaceDensity);
         if size(z, 2) ~= sz,
            error('Expected %d mass components but got %d.', ...
                  sz, size(z,2));
         end

         if numel(p) ~= size(z, 1),
            error('Got %d pressure values and %d mass values.', ...
                  numel(p), size(z,1));
         end

         i = 0;
         a = cell([1, max(1, nargout)]);
         for f = 1 : numel(self.fluids),
            sz = numel(self.fluids{f}.surfaceDensity);
            ix = i + 1 : i + sz;

            [a{:}] = self.fluids{f}.pvt(p, z(:,ix));

            for k = 1 : numel(a), varargout{k}(:,ix) = a{k}; end

            i = i + sz;
         end
      end
   end

   methods (Access=protected)
      function value = get_property(self, prop)
         %Retrieve value of given property for all constituent fluids.
         %
         % SYNOPSIS:
         %   value = self.get_property(prop);
         %
         % PARAMETERS:
         %   prop - Name of specific property for constituent fluids.
         %
         % RETURNS:
         %   value - Cell array of numerical values of this property.  One
         %           element for each constituent fluid component.

         value = cell([1, numel(self.fluids)]);
         for f = 1 : numel(value),
            value{f} = self.fluids{f}.(prop);
         end
      end
   end

   methods
      function value = get.surfaceDensity(self)
         value = self.get_property('surfaceDensity');
         value = [value{:}];
      end

      function value = get.info(self)
         value = strvcat(self.get_property('info'));  %#ok
      end

      function value = get.viscosity(self)
         value = self.get_property('viscosity');
         value = [value{:}];
      end

      function value = get.miscible(self)
         value = self.get_property('miscible');
         value = [value{:}];
      end
   end
end
