function G = cartGrid(celldim, varargin)
%Construct 2d or 3d Cartesian grid in physical space.
%
% SYNOPSIS:
%   G = cartGrid(celldim);
%   G = cartGrid(celldim, physdim);
%
% PARAMETERS:
%   celldim - Vector, length 2 or 3, specifying number of
%             cells in each coordinate direction.
%   physdim - Vector, length numel(celldim), of physical size in units of 
%             meters of the computational domain.
%             OPTIONAL.  Default value == celldim
%             (i.e., each cell has physical dimension 1-by-1-by-1 m).
%
% RETURNS:
%   G - Grid structure mostly as detailed in grid_structure, though lacking
%       the fields
%         - G.cells.volumes
%         - G.cells.centroids
%
%         - G.faces.areas
%         - G.faces.normals
%         - G.faces.centroids
%
%       These fields may be computed using the function computeGeometry.
%
%       There is, however, an additional field not described in
%       grid_structure:
%
%         - cartDims -- A length 3 vector giving number of cells in each
%                       coordinate direction.  In other words
%
%                               cartDims == celldim .
%
%       G.cellFaces(:,2) contains integers 1-6 corresponding to directions
%       W, E, S, N, T, B respectively.
%
% EXAMPLE:
%   % Make a 10-by-5-by-2 grid on the unit cube.
%      nx = 10; ny = 5; nz = 2;
%      G = cartGrid([nx, ny, nz], [1 1 1]);
%
%   % Plot the grid in 3D-view.
%      f = plotGrid(G); view(3);
%
% SEE ALSO:
%   grid_structure, tensorGrid, computeGeometry.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-12 12:12:22 +0200 (ma, 12 okt 2009) $
% $Revision: 2976 $

error(nargchk(1, 2, nargin, 'struct'));

if nargin < 2,
   physdim = celldim;
else
   physdim = varargin{1};
end

dx = repmat(physdim(1) ./ celldim(1), [celldim(1), 1]);
dy = repmat(physdim(2) ./ celldim(2), [celldim(2), 1]);
if numel(celldim) == 3,
   dz = repmat(physdim(3) ./ celldim(3), [celldim(3), 1]);
   G  = tensorGrid(dx, dy, dz);
else
   G  = tensorGrid(dx, dy);   
end

