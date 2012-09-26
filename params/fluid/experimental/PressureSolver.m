classdef PressureSolver < handle
    %
    methods(Abstract=true, Access=public)
       %
        [resSol, wellSol, S] = solve(self, resSol, wellSol, g, fluid)          
    end
end
