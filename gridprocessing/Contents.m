% Construct and manipulate reservoir simulation grid data structure.
%
% Files
%   boundaryFaceIndices - Retrieve face indices belonging to subset of global outer faces.
%   boundaryFaces       - Extract boundary faces from set of grid cells.
%   buildCornerPtNodes  - Construct physical nodal coordinates for CP grid.
%   cart2active         - Compute active cell numbers from linear Cartesian index.
%   cartGrid            - Construct Cartesian grid in physical space.
%   cellNodes           - Extract local-to-global vertex numbering for grid cells.
%   computeGeometry     - Compute geometry of grid.
%   computeGeometryS    - Compute geometry of grid using single precision.
%   deactivateZeroPoro  - Mark cells of zero porosity as inactive.
%   grid_structure      - Grid structure used in MATLAB Reservoir Simulation Toolbox.
%   processGRDECL       - Compute grid topology and geometry from pillar grid description.
%   removeCells         - Remove cells from grid and renumber cells, faces and nodes.
%   tensorGrid          - Construct Cartesian grid with variable physical cell sizes (dx,dy,dz).

%{
#COPYRIGHT#
%}
