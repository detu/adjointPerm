function varargout = outlineCoarseGrid(G, p, varargin)
%Impose outline of coarse grid on existing grid plot.
%
% SYNOPSIS:
%       outlineCoarseGrid(G, p)
%       outlineCoarseGrid(G, p, 'pn1', pv1, ...)
%   h = outlineCoarseGrid(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   p       - Coarse grid partition vector as defined by (e.g) partitionUI
%             and processPartition.
%
%   'pn'/pv - List of PATCH property specifcations.
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.
%
% EXAMPLE:
%   G  = cartGrid([8, 8, 2]);
%   p  = partitionUI(G, [2, 2, 1]);
%   % plot fine grid:
%   plotGrid(G, 'faceColor', 'none'); view(3);
%   % outline coarse grid on fine grid:
%   outlineCoarseGrid(G, p);
%
% SEE ALSO:
%   plotFaces, plotGrid, patch.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-14 09:35:58 +0200 (on, 14 okt 2009) $
% $Revision: 2997 $

   dim = size(G.nodes.coords,2);
   assert (sum(dim == [2, 3]) == 1);

   if dim == 2,
      % Include outer boundary in 2D.
      mask = true([G.faces.num, 1]);
   else
      mask = all(G.faces.neighbors > 0, 2);
   end

   p = [0; p];
   N = p(G.faces.neighbors + 1);
   N = find((N(:,1) ~= N(:,2)) & mask);

   if dim == 2,
      h = plotEdges(G, N, [0.85, 0.15, 0.85]);
   else
      h = plotFaces(G, N, [0.85, 0.15, 0.85], varargin{:});
   end

   if nargout > 0, varargout{1} = h; end
end

%--------------------------------------------------------------------------

function h = plotEdges(G, e, c)
   assert (all(diff([G.faces.nodePos(  e  ), ...
                     G.faces.nodePos(e + 1)], [], 2) == 2));

   fn = G.faceNodes(mcolon(G.faces.nodePos(  e  ), ...
                           G.faces.nodePos(e + 1) - 1));
   fc = G.nodes.coords(fn,:);

   xd = reshape(fc(:,1), 2, []);
   yd = reshape(fc(:,2), 2, []);
   h  = line(xd, yd, 'Color', c, 'LineWidth', 1);
end
