function varargout = plotCellData(G, data, varargin)
%Plot exterior grid faces, coloured by given data, to current axes.
%
% SYNOPSIS:
%       plotCellData(G, data)
%       plotCellData(G, data, 'pn1', pv1, ...)
%       plotCellData(G, data, cells)
%       plotCellData(G, data, cells, 'pn1', pv1, ...)
%   h = plotCellData(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   data    - Scalar cell data with which to colour the grid.  One scalar
%             value for each cell in the grid.  If a cell subset is
%             specified in terms of the 'cell' argument, 'data' must either
%             contain one scalar value for each cell in the model or one
%             scalar value for each cell in this subset.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function 'plotCellData' will be restricted to the
%             subset of cells from 'G' represented by 'cells'.
%
%             If empty or unspecified, function 'plotCellData' will behave
%             as if the user defined
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
%   Function 'plotCellData' is implemented directly in terms of the
%   low-level function PATCH.  If a separate axes is needed for the
%   graphical output, callers should employ function newplot prior to
%   calling 'plotCellData'.
%
% EXAMPLE:
%   Given a grid 'G' and a reservoir solution structure 'resSol' returned
%   from, e.g., function 'solveIncompFlow', plot the cell pressure in bar:
%
%      figure, plotCellData(G, convertTo(resSol.cellPressure, barsa()));
%
% SEE ALSO:
%   plotFaces, boundaryFaces, patch, newplot.

%{
#COPYRIGHT#
%}

% $Id: plotCellData.m 2338 2009-06-05 17:19:30Z bska $

if numel(G) > 1,
   error(msgid('Grid:MultiComponent'), ...
         'Cannot plot more than one grid at a time.');
end

% Default to providing graphical output from all cells in the grid model.
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
assert (numel(data) == G.cells.num || numel(data) == numel(cells));

if size(G.nodes.coords, 2) == 3,
   [f, c] = boundaryFaces(G, cells);
else
   % For 2D grids, the faces to plot are the actual individual grid cells.
   [f, c] = deal(cells);
end

if numel(data) < G.cells.num,
   renum        = zeros([G.cells.num, 1]);
   renum(cells) = 1 : numel(cells);
   c            = renum(c);

   assert (all(c > 0) && all(c <= numel(data)));
end
h = plotFaces(G, f, data(c), 'EdgeColor', 'none', varargin{:});

if nargout > 0, varargout{1} = h; end
