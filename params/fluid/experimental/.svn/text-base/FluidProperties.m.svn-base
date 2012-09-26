classdef FluidProperties
    properties(Access=public)
        fluid
    end
    methods(Access=public)
        function self = FluidProperties(f)
            self.fluid = f;
        end

        function m = mobility(self, s, p)
            [d,c,v] = self.fluid.getPressureDependentData(p);
            kr      = self.fluid.getSatDependentData(s);

            m = kr./v;
        end

        function lt = totalMobility(self, s, p)
            lt = sum(self.mobility(s, p), 2);
        end
    end
end