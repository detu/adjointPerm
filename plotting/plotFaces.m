function varargout = plotFaces(G, faces, varargin)
%Plot selection of coloured grid faces to current axes (reversed Z axis).
%
% SYNOPSIS:
%       plotFaces(G, faces)
%       plotFaces(G, faces, 'pn1', pv1, ...)
%       plotFaces(G, faces, colour)
%       plotFaces(G, faces, colour, 'pn1', pv1, ...)
%   h = plotFaces(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   faces   - Vector of face indices.  The graphical output of 'plotFaces'
%             will be restricted to the subset of grid faces from 'G'
%             represented by 'faces'.
%
%   colour  - Colour data specification.  Either a MATLAB 'ColorSpec'
%             (i.e., an RGB triplet (1-by-3 row vector) or a short or long
%             colour name such as 'r' or 'cyan'), or a PATCH
%             'FaceVertexCData' table suiteable for either indexed or
%             'true-colour' face colouring.  This data *MUST* be an m-by-1
%             column vector or an m-by-3 matrix.  We assume the following
%             conventions for the size of the colour data:
%
%                - ANY(SIZE(colour,1) == [1, NUMEL(faces)])
%                  One (constant) indexed colour for each face in 'faces'.
%                  This option supports 'flat' face shading only.  If
%                  SIZE(colour,1) == 1, then the same colour is used for
%                  all faces in 'faces'.
%
%                - SIZE(colour,1) == G.nodes.num
%                  One (constant) indexed colour for each node in 'faces'.
%                  This option must be chosen in order to support
%                  interpolated face shading.
%
%             OPTIONAL.  Default value: colour = 'y' (shading flat).
%
%   'pn'/pv - List of other property name/value pairs.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h - Handle to resulting PATCH object.  The patch object is added to the
%       current AXES object.
%
% NOTES:
%   Function 'plotFaces' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotFaces'.
%
% EXAMPLE:
%   % Plot grid with boundary faces on left side in red colour:
%   G     = cartGrid([5, 5, 2]);
%   faces = boundaryFaceIndices(G, 'LEFT', 1:5, 1:2, []);
%   plotGrid (G, 'faceColor', 'none'); view(3)
%   plotFaces(G, faces, 'r');
%
% SEE ALSO:
%   plotCellData, plotGrid, newplot, patch, shading.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-13 16:55:17 +0200 (ti, 13 okt 2009) $
% $Revision: 2992 $

   if isempty(faces), return, end

   if ~isnumeric(faces),
      error(msgid('FaceList:NonNumeric'), ...
            'Face list ''faces'' is not numeric.')
   end

   [plotOutline, varargin] = do_outline_p(varargin{:});

   if mod(numel(varargin), 2) == 0,
      colour   = 'yellow';
   else
      colour   = varargin{1};
      varargin = { varargin{2 : end} };
   end

   % Extract face topology for subset of grid faces and the coordinates for
   % the actual vertices ('verts') present in this topology.
   %
   dim = size(G.nodes.coords, 2);   assert (any(dim == [2, 3]));
   face_topo  = str2func(['get_face_topo_', int2str(dim), 'd']);
   [f, verts] = face_topo     (G, faces);
   v          = G.nodes.coords(verts, :);

   % Massage colour data into form suitable for 'FaceVertexCData' property.
   %
   % From here, we assume that 'colour' is an m-by-1 (for indexed
   % colouring) or m-by-3 (for 'true colour' colouring) array.
   % Furthermore, 'm' is assumed to be either 1 (meaning all 'faces' should
   % be coloured using a single colour), NUMEL(faces) when the individual
   % faces are coloured separately, or G.nodes.num for interpolated face
   % colouring.
   %
   if ischar(colour), colour = get_rgb(colour); end

   if any(size(colour,1) == [1, numel(faces)]),
      fc     = 'flat';
   elseif size(colour,1) == G.nodes.num,
      colour = colour(verts,:);
      fc     = 'interp';
   end

   % Build final patch for graphical output (Note: added to GCA).
   h = patch('Faces'          , f      , 'Vertices' , v , ...
             'FaceVertexCData', colour , 'FaceColor', fc, varargin{:});

   set(get(h, 'Parent'), 'ZDir', 'reverse')

   if plotOutline,
      pts = findFaceOutline(G, faces);
      do_hold = ishold();
      hold on, plot3(pts(:,1), pts(:,2), pts(:,3), 'k');
      if ~do_hold, hold off, end
   end

   if nargout > 0, varargout{1} = h; end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [f, present] = get_face_topo_2d(G, cells)  %#ok
   eIX = G.cells.facePos;
   nn  = double(diff([G.cells.facePos(cells), ...
                      G.cells.facePos(cells + 1)], [], 2));
   cellNodes = getCellNodes(G);
   cn  = double(cellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));

   m   = numel(cells);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %
   present           = false([G.nodes.num, 1]);
   present(cn)       = true;

   node_num          = zeros([G.nodes.num, 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(cn);

   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
end

%--------------------------------------------------------------------------

function [f, present] = get_face_topo_3d(G, faces)  %#ok
   eIX = G.faces.nodePos;
   nn  = double(diff([G.faces.nodePos(faces), ....
                      G.faces.nodePos(faces + 1)], [], 2));
   fn  = double(G.faceNodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

   m   = numel(faces);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %
   present           = false([G.nodes.num, 1]);
   present(fn)       = true;

   node_num          = zeros([G.nodes.num, 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(fn);

   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
end

%--------------------------------------------------------------------------

function [plotOutline, varargin] = do_outline_p(varargin)
   % Does caller of 'plotFaces' request an outline plot?

   plotOutline = false;
   if numel(varargin) > 0,
      i = 1 + isnumeric(varargin{1});
      v = find(strcmpi({ varargin{i : 2 : end} }, 'outline'));
      if ~isempty(v),
         % Find argument following last 'outline'
         plotOutline = varargin{i + 2*v(end) - 1};

         % Remove pairs of 'outline'/value from varargin.
         varargin([i+2*v-2, i+2*v-1]) = [];
      end
   end
end

%--------------------------------------------------------------------------

function rgb = get_rgb(colour)
   switch lower(colour),
      case {'y', 'yellow' }, rgb = [1, 1, 0];
      case {'m', 'magenta'}, rgb = [1, 0, 1];
      case {'c', 'cyan'   }, rgb = [0, 1, 1];
      case {'r', 'red'    }, rgb = [1, 0, 0];
      case {'g', 'green'  }, rgb = [0, 1, 0];
      case {'b', 'blue'   }, rgb = [0, 0, 1];
      case {'w', 'white'  }, rgb = [1, 1, 1];
      case {'k', 'black'  }, rgb = [0, 0, 0];
      otherwise            , rgb = [0, 0, 1]; % Unknown colour -> 'blue'.
   end
end

%--------------------------------------------------------------------------
% For 2d only.
function cn = getCellNodes(G)
   % Construct n x 2 table of cell edges with edges oriented the same
   % direction around the cell boundary.
   cellNo    = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
   edges     = reshape(G.faceNodes, 2, [])';
   cellEdges = edges(G.cellFaces(:,1),:);
   ind       = G.faces.neighbors(G.cellFaces(:,1), 1) ~= cellNo;
   cellEdges(ind, :) = cellEdges(ind, [2,1]);

   % Sort edges in each cell:
   for c = 1 : G.cells.num,
      ind = G.cells.facePos(c) : G.cells.facePos(c + 1) - 1;
      cellEdges(ind, :) = sortEdges(cellEdges(ind,:));
   end
   cn = reshape(cellEdges(:,1), 1, [])';
end

%--------------------------------------------------------------------------
% For 2d only.
function edges = sortEdges(edges)
   % Assume edges vectors are oriented in the same direction around cell.
   % Sort edges such that they are back-to-back.
   % Then cellNodes are edges(:,1).

   for i = 1 : size(edges, 1) - 1,
      for j = i + 1 : size(edges,1),
         if edges(i,2) == edges(j,1), break; end
      end

      % Swap edges i+1 and j
      tmp = edges(i+1,:);
      edges(i+1, :) = edges(j,:);
      edges(j,   :) = tmp;
   end
end

%--------------------------------------------------------------------------

function pts = findFaceOutline(g, faces)
   assert (size(g.nodes.coords, 2) == 3);
   % Find points on border of collection of grid faces.
   if numel(faces)==0, pts = zeros(0,3); return; end

   cellNo = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
   sgn    = 2*(cellNo == g.faces.neighbors(g.cellFaces(:,1), 1)) - 1;
   faceNo = rldecode(1:g.faces.num, diff(g.faces.nodePos), 2) .';

   fi     = false(g.faces.num, 1); fi(faces) = true;
   fn     = g.faceNodes(fi(faceNo));

   nodeflag     = false([g.nodes.num, 1]);
   nodeflag(fn) = true;

   fe = faceEdges(g);
   fe(sgn(faceNo)<0, :) = fe(sgn(faceNo)<0, [2,1]);

   fe = fe (fi(faceNo),:);
   nodes = double(fe(all(nodeflag(fe), 2), :));

   if numel(nodes) > 0,
      % Remove edges which appear more than once.  These edges are
      % considered internal in the collection of faces.  The remaining
      % edges are on the outer boundary.
      [nodes, n] = rlencode(sortrows(sort(nodes, 2)));
      nodes(n>1,:) = [];
   end

   pts = nan(size(nodes, 1)*3, 3);
   pts (1:3:end,:) = g.nodes.coords(nodes(:,1),:);
   pts (2:3:end,:) = g.nodes.coords(nodes(:,2),:);
end

%--------------------------------------------------------------------------

function fe = faceEdges(g)
   fe = [g.faceNodes, g.faceNodes([2:end,1])];
   fe(g.faces.nodePos(2:end)-1,2) = fe(g.faces.nodePos(1:end-1),1);
end
