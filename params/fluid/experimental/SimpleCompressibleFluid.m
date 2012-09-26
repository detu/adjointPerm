classdef SimpleCompressibleFluid < Fluid
    % Simple compressible two-phase fluid model for a water-oil system.
    % The water and oil viscosities are 0.1 and 1.0 cP, and the densities
    % are 1014 and 780 Kg/m3 at 100 bar.  The oil compressibility is 0.006.
    %
    % The relative permeabilities are simple quadratic functions.

    properties(Access=protected)
       reference_pressure    = 100.0;
       reference_density     = [1014, 780.0];

       compressibility       = [0.0, 0.006];
       viscosity             = [0.1, 1.0  ];

    end


    methods(Access=public)
        function self = SimpleCompressibleFluid()
        end

        function kr = getSatDependentData(self, s)
            assert(size(s,2)==1,  'Saturation must be a column vector')
            kr  = [s.^2, (1 - s).^2];
        end


        function [d, c, v] = getPressureDependentData(self, p)
            n   = numel(p);
            t   = bsxfun(@times, p(:)-self.reference_pressure, self.compressibility);

            d = bsxfun(@times, self.reference_density, exp(t));
            c = repmat(self.compressibility, [n,1]);
            v = repmat(self.viscosity,       [n,1]);
        end
    end
end