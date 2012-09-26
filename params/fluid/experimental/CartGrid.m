classdef CartGrid < TensorGrid
    methods(Access=public)
        function self = CartGrid(celldim, varargin)
            % CARTGRID -- Construct cartesian grid in physical space
            %
            % SYNOPSIS:
            %   G = cartGrid(celldim);
            %   G = cartGrid(celldim, physdim);
            %
            % PARAMETERS:
            %   celldim - Vector, length 3, specifying number of
            %             cells in each coordinate direction.
            %   physdim - Vector, length 3, of physical dimensions.
            %             OPTIONAL.  Default value == celldim
            %             (i.e. each cell has physical dimension 1-by-1-by-1).
            %
            % RETURNS:
            %   G - Grid structure mostly as detailed in grid_structure, though lacking the
            %       fields
            %         - G.cells.volumes
            %         - G.cells.centroids
            %         - G.cells.ijkMap
            %
            %         - G.faces.areas
            %         - G.faces.normals
            %         - G.faces.centroids
            %
            %       These fields may be computed using the function computeGeometry.
            %
            %       There is, however, an additional field not described in grid_structure:
            %
            %         - cartDims -- A length 3 vector giving number of cells in each
            %                       coordinate direction.  In other words
            %
            %                               cartDims == celldim .
            %
            % SEE ALSO:
            %   grid_structure, tensorGrid, computeGeometry.

            % $Id: cartGrid.m 153 2008-04-29 13:19:15Z jrn $

            error(nargchk(1, 2, nargin, 'struct'));

            if nargin < 2,
                physdim = celldim;
            else
                physdim = varargin{1};
            end

            dx = repmat(physdim(1) ./ celldim(1), [celldim(1), 1]);
            dy = repmat(physdim(2) ./ celldim(2), [celldim(2), 1]);
            dz = repmat(physdim(3) ./ celldim(3), [celldim(3), 1]);

            self  = self@TensorGrid(dx, dy, dz);
        end
    end
end