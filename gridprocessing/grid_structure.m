%Grid structure used in MATLAB Reservoir Simulation Toolbox.
%
% SYNOPSIS:
%   1) Construct Cartesian grid.
%        G = cartGrid(...);
%        G = computeGeometry(G);
%
%   2) Read corner point grid specification from file.
%        grdecl = readGRDECL(...);
%        G      = processGRDECL(grdecl);
%        G      = computeGeometry(G);
%
% RETURNS:
%   G - Grid structure.  A master structure having the following fields:
%         - cells      --  A structure specifying properties for each
%                          individual cell in the grid.  See CELLS below
%                          for details.
%
%         - faces      --  A structure specifying properties for each
%                          individual (global) face in the grid.  See FACES
%                          below for details.
%
%         - nodes      --  A structure specifying properties for each
%                          individual node (vertex) in the grid.  See NODES
%                          below for details.
%
%         - cellFaces  --  A (G.cells.facePos(end)-1)-by-3 array of global
%                          faces connected to a given cell.  Specifically,
%                          if cellFaces(i,1)==j, then global face
%                          cellFaces(i,2) is connected to global cell `j'
%                          and cellFaces(i,3) contains a tag by which local
%                          face direction (East, West, South, North,
%                          Bottom, Top) may be distinguished.  To conserve
%                          memory, only the final two columns are stored.
%                          The first column may be re-constructed using the
%                          statement
%                              rldecode(1 : G.cells.num, ...
%                                       diff(G.cells.facePos), 2) .'
%
%         - faceNodes  --  A (G.faces.nodePos(end)-1)-by-2 array of vertices
%                          in the grid.  Specifically, if
%                          faceNodes(i,1)==j, then local vertex `i' is part
%                          of global face number `j' and corresponds to
%                          global vertex faceNodes(i,2).  To conserve
%                          memory, only the last column is stored.  The
%                          first column can be constructed using the
%                          statement
%                              rldecode(1 : G.faces.num, ...
%                                       diff(G.faces.nodePos), 2) .'
%
%  CELLS - Cell structure G.cells:
%    - num       --  Number of cells in global grid.
%
%    - numFaces  --  Number of half-faces for each cell in global grid
%                    (== REPMAT(6, [G.cells.num, 1]) for Cartesian geom.).
%
%                    The face information for cell `i' is found in the
%                    submatrix cellFaces(IX(i)+1 : IX(i+1), :) provided
%
%                          IX = CUMSUM([0; double(numFaces)]) .
%
%                    Note: This field is deprecated and will be removed.
%
%    - facePos   --  Indirection map of size [num+1,1] into the 'cellFaces'
%                    array.  Specifically, the face information of cell 'i'
%                    is found in the submatrix
%
%                          cellFaces(facePos(i) : facePos(i+1)-1, :)
%
%                    The number of faces of each cell may be computed using
%                    the statement DIFF(facePos).
%
%    - indexMap  --  Maps internal to external grid cells (i.e., active
%                    cell numbers to global cell numbers).  In the case of
%                    Cartesian grids, indexMap == (1 : G.cells.num)' .
%
%                    For logically Cartesian grids of Cartesian dimension
%                    'dims', e.g., corner point grids, a map of cell
%                    numbers to logical indices may be constructed using
%                    the statements
%
%                       [ijk{1:3}] = ind2sub(dims, G.cells.indexMap(:));
%                       ijk        = [ijk{:}];
%
%                    Then ijk(i,:) is the global (I,J,K) index of cell 'i'.
%
%    - volumes   --  A G.cells.num-by-1 array of cell volumes.
%
%    - centroids --  A G.cells.num-by-3 array of cell centroids.
%
%  FACES - Face structure G.faces:
%    - num       --  Number of global faces in grid.
%
%    - numNodes  --  Number of half-edges for each face in global grid
%                    (== REPMAT(4, [G.faces.num, 1])) for Cartesian geom.).
%
%                    The node information for global face `i' is found in
%                    the submatrix faceNodes(IX(i)+1 : IX(i+1), :) provided
%
%                          IX = CUMSUM([0; double(numNodes)]) .
%
%                    Note: This field is deprecated and will be removed.
%
%    - nodePos   --  Indirection map of size [num+1,1] into the 'faceNodes'
%                    array.  Specifically, the node information of face 'i'
%                    is found in the submatrix
%
%                          faceNodes(nodePos(i) : nodePos(i+1)-1, :)
%
%                    The number of nodes of each face may be computed using
%                    the statement DIFF(nodePos).
%
%    - neighbors --  A G.faces.num-by-2 array of neighbouring information.
%                    Global face `i' is shared by global cells neigbors(i,1)
%                    and neighbors(i,2).  One of neighbors(i,1) or
%                    neighbors(i,2), but not both, may be zero, meaning that
%                    face `i' is an external face shared only by a single cell
%                    (the non-zero one).
%
%    - areas     --  A G.faces.num-by-1 array of face areas.
%
%    - normals   --  A G.faces.num-by-3 array of *AREA WEIGHTED*, directed
%                    face normals.  The normal on face `i' points from cell
%                    G.faces.neighbors(i,1) to cell G.faces.neighbors(i,2).
%
%    - centroids --  A G.faces.num-by-3 array of face centroids.
%
%  NODES - Node (vertex) structure G.nodes:
%    - num       --  Number of global nodes (vertices) in grid.
%
%    - coords    --  A G.nodes.num-by-3 array of physical nodal coordinates.
%                    Global node `i' is at physical coordinate [coords(i,1:3)]
%
% REMARKS:
%   The grid is constructed according to a right-handed coordinate system
%   where the Z coordinate is interpreted as depth.  Consequently, plotting
%   routines such as plotGrid display the grid with a reversed Z axis.
%
% SEE ALSO:
%   In case the description scrolled off screen too quickly, you may access
%   the information at your own pace using the command
%
%               more on, help grid_structure, more off

%{
#COPYRIGHT#
%}

% $Date: 2009-09-28 16:21:53 +0200 (ma, 28 sep 2009) $
% $Revision: 2870 $
