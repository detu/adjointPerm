classdef Fluid
    % Base for fluid objects implementing common functionality
    % Cannot be instantiated.


    % Public interfaces
    methods(Access=public, Abstract=true)
        [d, c, v] = getPressureDependentData (self, p)
        kr        = getSatDependentData      (self, s)
    end
end
