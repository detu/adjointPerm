function [f, varargout] = boundaryFaces(g, varargin)
%Extract boundary faces from set of grid cells.
%
% SYNOPSIS:
%    f     = boundaryFaces(G)
%    f     = boundaryFaces(G, cells)
%   [f, c] = boundaryFaces(...)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   cells - Non-empty subset of cells from which to extract boundary faces.
%           OPTIONAL.  Default value: cells = 1 : G.cells.num, meaning all
%           external faces for all grid cells will be extracted.  This
%           amounts to extracting the entire boundary of 'G'.
%
% RETURNS:
%   f - List of faces bounding the sub domain given by 'cells'.
%
%   c - List of specific grid cells connected to the individual faces in
%       'f'.  This may be useful for plotting cell data (e.g., the cell
%       pressure) on the sub domain faces by means of function 'plotFaces'.
%       OPTIONAL.  Only returned if specifically requested.
%
% EXAMPLE:
%   G    = cartGrid([40, 40, 5]);
%   rock = <some rock data structure for G>;
%
%   % 1) Plot (external) geometry of 'G'.
%   f  = boundaryFaces(G);
%   hg = plotFaces(G, f);
%
%   % 2) Plot horizontal permeability along diagonal of reservoir
%   [f, c] = boundaryFaces(G, 1 : G.cartDims(1) + 1 : G.cells.num);
%   hd     = plotFaces(G, f, log10(rock.perm(c,1)), 'FaceAlpha', 0.625);
%
% SEE ALSO:
%   plotFaces.

%{
#COPYRIGHT#
%}

% $Id: boundaryFaces.m 1949 2009-03-31 10:16:08Z bska $

   sub = (1 : g.cells.num) .';
   if (nargin > 1) && isnumeric(varargin{1}),
      sub = varargin{1};
   end

   present          = false([g.cells.num + 1, 1]);
   present(sub + 1) = true;

   % Extract boundary faces:
   %   Those faces for which one of the connecting cells is 'present'
   %   (i.e., within the sub domain 'sub') and the other is ~present
   %   (outside the sub domain).  This means that
   %
   %       SUM(present(g.faces.neighbors), 2) == 1
   %
   %   for the external faces.  Add one to account for the outside being
   %   'cell 0'.
   %
   f = find(sum(double(present(g.faces.neighbors + 1)), 2) == 1);

   if nargout > 1,
      % User requested list of cells connected to the faces in 'f'.

      % We construct a column vector 'c' such that c(i) is the cell within
      % 'sub' which connects to the face f(i) for all i=1:numel(f).
      %
      % Algorithm:
      %   1) Construct reduced connection matrix 'n', one row for each face
      %      in 'f', such that one of the cells n(i,1) or n(i,2) is the
      %      required cell c(i).
      %
      %   2) Mask reduced connection matrix by cell presence.  Result is a
      %      matrix whose entries are either zero (where the original
      %      entries of 'n' contained a cell outside 'sub') or the original
      %      cell number (if inside 'sub').
      %
      %   3) Sum the rows of this matrix to obtain the corresponding cell
      %      number regardless of original column number.

      n = double(g.faces.neighbors(f,:));
      c = sum(n .* double(present(n + 1)), 2);

      assert (all(c > 0));

      varargout{1} = c;
   end
end
