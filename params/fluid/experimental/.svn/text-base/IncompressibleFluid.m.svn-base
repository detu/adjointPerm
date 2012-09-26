classdef IncompressibleFluid < Fluid
    % Base for incompressible fluids.  Cannot be instantiated.

    properties(Access=protected)
       viscosity
       compressibility
       density
    end
    
    % Member functions that need to be implemented in subclasses
    methods(Access=public)
        

        function [d, c, v] = getPressureDependentData(self, p)
            v = repmat(self.viscosity,       [numel(p), 1]);
            c = repmat(self.compressibility, [numel(p), 1]);
            d = repmat(self.density,         [numel(p), 1]);
        end
    end
end

