classdef SimpleIncompressibleFluid < IncompressibleFluid
    % Simple incompressible two-phase fluid model for a water-oil system.
    % The water and oil viscosities are 0.1 and 1.0 cP, and the densities
    % are 1014 and 780 Kg/m3.
    %
    % The relative permeabilities are simple quadratic functions.
    
    methods(Access=public)
        function self = SimpleIncompressibleFluid()
            self.density         = [1014.0, 780.0];
            self.viscosity       = [   0.1,   1.0];
            self.compressibility = [   0.0,   0.0];
        end

    end
    methods(Access=public)
        function kr = getSatDependentData(self, s)
            assert(size(s,2)==1,  'Saturation must be a column vector')
            kr  = [s.^2, (1 - s).^2];
        end
    end
end