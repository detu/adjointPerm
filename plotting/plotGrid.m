function varargout = plotGrid(G, varargin)
%Plot exterior grid faces to current axes (reversed Z axis).
%
% SYNOPSIS:
%       plotGrid(G)
%       plotGrid(G, 'pn1', pv1, ...)
%       plotGrid(G, cells)
%       plotGrid(G, cells, 'pn1', pv1, ...)
%   h = plotGrid(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function 'plotGrid' will be restricted to the
%             subset of cells from 'G' represented by 'cells'.
%
%             If empty or unspecified, function 'plotGrid' will behave as
%             if the caller defined
%
%                 cells = 1 : G.cells.num
%
%             meaning graphical output will be produced for all cells in
%             the grid model 'G'.
%
%   'pn'/pv - List of property names/property values.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current AXES object (GCA).
%        OPTIONAL.  Only returned if specifically requested.
%
% NOTES:
%   Function 'plotGrid' is implemented directly in terms of the low-level
%   function PATCH.  If a separate axes is needed for the graphical output,
%   callers should employ function newplot prior to calling 'plotGrid'.
%
% EXAMPLE:
%   G = cartGrid([10, 10, 5]);
%
%   % 1) Plot grid with yellow colour on faces (default):
%   figure, plotGrid(G, 'EdgeAlpha', 0.1); view(3)
%
%   % 2) Plot grid with no colour on faces (transparent faces):
%   figure, plotGrid(G, 'FaceColor', 'none'); view(3)
%
% SEE ALSO:
%   plotCellData, plotFaces, patch, newplot.

%{
#COPYRIGHT#
%}

% $Id: plotGrid.m 2338 2009-06-05 17:19:30Z bska $

if numel(G) > 1,
   error(msgid('Grid:MultiComponent'), ...
         'Cannot plot more than one grid at a time.');
end

% Default to showing the entire grid model.
%
cells = (1 : G.cells.num) .';

if mod(numel(varargin), 2) == 1,
   % Caller requested graphical output from a particular subset of the grid
   % cells.  Honour request, but treat an empty 'cells' argument as if the
   % caller requested the default behaviour (i.e., cells = 1:G.cells.num).
   %
   if isnumeric(varargin{1}) && ~isempty(varargin{1}),
      cells = varargin{1};
   end

   % Strip 'cells' argument off of remaining input argument list.
   %
   varargin = { varargin{2 : end} };
end

% Assert that remaining arguments at least appear to be 'name'/value pairs
% intended for PATCH.
%
assert (all(cellfun(@ischar, { varargin{1 : 2 : end} })));

if size(G.nodes.coords, 2) == 3,
   f = boundaryFaces(G, cells);
else
   % For 2D grids, the faces to plot are the actual individual grid cells.
   f = cells;
end

h = plotFaces(G, f, 'EdgeColor', 'k', 'FaceColor', 'y', varargin{:});

if nargout > 0, varargout{1} = h; end
