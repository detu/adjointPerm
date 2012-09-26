classdef SinglePhaseIncompressibleFluid < IncompressibleFluid
    % Simple incompressible water model.
    % The water viscosity is 0.1 cP, and the density is 1014 Kg/m3.
    %
    % The relative permeability is 1.0.
    
    methods(Access=public)
        function self = SinglePhaseIncompressibleFluid()
            self.density         = [1014.0];
            self.viscosity       = [   0.1];
            self.compressibility = [   0.0];
        end

    end
    methods(Access=public)
        function kr = getSatDependentData(self, s)
            assert(size(s,2)==1,  'Saturation must be a column vector')
            kr  = ones(size(s));
        end
    end
end