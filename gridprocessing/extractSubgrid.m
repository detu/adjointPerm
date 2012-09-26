function H = extractSubgrid(G, c)
%Construct valid grid definition from subset of existing grid cells.
%
% SYNOPSIS:
%   subG = extractSubgrid(G, cells)
%
% PARAMETERS:
%   G     - Valid grid definition.
%
%   cells - Subset of existing grid cells for which the subgrid data
%           structure should be constructed.
%
% RETURNS:
%   subG  - Resulting subgrid.
%
% EXAMPLE:
%   G    = cartGrid([3, 5, 7]);     [I, J, K] = ndgrid(2, 2:4, 3:5);
%   subG = extractSubgrid(G, sub2ind(G.cartDims, I(:), J(:), K(:)));
%   plotGrid(subG, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.25)
%   view(-35,20), camlight
%
% SEE ALSO:
%   grid_structure.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-23 18:33:51 +0200 (ti, 23 jun 2009) $
% $Revision: 2398 $

   % We don't support multi-component grids...
   if numel(G) > 1,
      error(msgid('Grid:MultiComponent'), ...
            'Cannot extract sub-grid from more than one grid at a time.');
   end

   if ~isfield(G.cells, 'facePos'),
      error(msgid('Grid:Cells:NoFacePos'), ...
           ['Grid must support ''facePos'' field in order to ', ...
            'extract sub-grid.']);
   end

   if ~isfield(G.faces, 'nodePos'),
      error(msgid('Grid:Faces:NoNodePos'), ...
           ['Grid must support ''nodePos'' field in order to ', ...
            'extract sub-grid.']);
   end

   ix = @(p,i) mcolon(double(p(i)), double(p(i+1)) - 1);

   % Sort the subset cells, c, in ascending order.  Add one for outside.
   cells      = false([G.cells.num + 1, 1]);
   cells(c+1) = true;

   % Extract faces connected to 'c'.
   faces = false([G.faces.num, 1]);
   faces(G.cellFaces(ix(G.cells.facePos, c))) = true;

   % Extract nodes connected to 'faces' (in cell subset).
   nodes = false([G.nodes.num, 1]);
   ff    = find(faces);
   nodes(G.faceNodes(ix(G.faces.nodePos, ff))) = true;

   H.nodes.coords    = G.nodes.coords(nodes, :);
   H.cells.numFaces  = G.cells.numFaces(cells(2:end));
   H.faces.numNodes  = G.faces.numNodes(faces);
   H.faces.neighbors = G.faces.neighbors(faces, :);

   pos = @(n) uint32(cumsum([1; double(reshape(n, [], 1))]));
   H.cells.facePos   = pos(H.cells.numFaces);
   H.faces.nodePos   = pos(H.faces.numNodes);

   H.cellFaces       = G.cellFaces(ix(G.cells.facePos, c ), :);
   H.faceNodes       = G.faceNodes(ix(G.faces.nodePos, ff)   );

   % Renumbering of grid entities (cells,faces,nodes) in subgrid.
   mc = zeros([G.cells.num+1, 1], 'uint32');   mc(cells) = 1:numel(c);
   fc = zeros([G.faces.num  , 1], 'uint32');   fc(faces) = 1:numel(ff);
   nc = zeros([G.nodes.num  , 1], 'uint32');   nc(nodes) = 1:sum(nodes);

   H.faces.neighbors = mc(H.faces.neighbors+1);
   H.cellFaces(:,1)  = fc(H.cellFaces(:,1));
   H.faceNodes       = nc(H.faceNodes);

   H.faces.num = numel(H.faces.numNodes);
   H.cells.num = numel(H.cells.numFaces);
   H.nodes.num = size(H.nodes.coords, 1);

   H.cells.indexMap = G.cells.indexMap(cells(2:end));

   H.cartDims = nan(size(G.cartDims));

   % Preserve 'computeGeometry' fields, if present.
   if isfield(G.cells, 'volumes'),
      H.cells.volumes   = G.cells.volumes  (cells(2:end)  );
   end
   if isfield(G.cells, 'centroids'),
      H.cells.centroids = G.cells.centroids(cells(2:end),:);
   end
   if isfield(G.faces, 'areas'),
      H.faces.areas     = G.faces.areas    (faces  );
   end
   if isfield(G.faces, 'normals'),
      H.faces.normals   = G.faces.normals  (faces,:);
   end
   if isfield(G.faces, 'centroids'),
      H.faces.centroids = G.faces.centroids(faces,:);
   end
end
