classdef Boundary < handle

    properties(Access=public)
       grid
       data  % n x 3 table of dirichlet/neumann marker, pressure/flux, saturation.
             % n=grid.faces.num. 1=Dirichlet, 2=Neumann.

    end

    methods(Access=public)
        function self = Boundary(grid)
            self.grid = grid;
            self.data = zeros(grid.faces.num, 4);

            % Homogeneous Neumann condition
            ind = prod(double(grid.faces.neighbors), 2) > 0;
            self.data(ind,1) = 2;  % Neumann condition.
        end

        function self = setDirichletCondition(self)
        end

        function self = setNeumannCondition(self)
        end

        function self = pside(self, side, I1, I2, pressure)
            ix = boundaryFaceIndices(self.grid, side, I1, I2);
            self.data(ix, 1) = 1;
            self.data(ix, 2) = pressure;
            self.data(ix, 3) = 1.0;
        end

        function self = fluxside(self)
        end


        function [faces, values, sat] = getDirichletBoundary(self)
            faces  = find(self.data(:,1) == 1);
            values = self.data(faces,2);
            sat    = self.data(faces,3);
        end

        function [faces, values, sat] = getNeumannBoundary(self)
            faces  = find(self.data(:,1) == 2);
            values = self.data(faces,2);
            sat    = self.data(faces,3);
        end

        function [faces, sat] = getBoundarySaturation(self)
            faces = find(self.data(:,1) > 0);
            sat   = self.data(ind,3);
        end
    end
end