function activeCells = cart2active(G, cartCells)
%Compute active cell numbers from linear Cartesian index.
%
% SYNOPSIS:
%   activeCells = cart2active(G, c)
%
% PARAMETERS:
%   G - Grid data structure as described by 'grid_structure'.
%   c - List of linear Cartesian cell indices.
%
% RETURNS:
%   activeCells - Active cell numbers corresponding to the individual
%                 Cartesian cell numbers in 'c'.
%
% NOTE:
%   This function provides the inverse mapping of the 'G.cells.indexMap'
%   field in the grid data structure.
%
% SEE ALSO:
%   grid_structure.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-28 19:17:34 +0200 (ma, 28 sep 2009) $
% $Revision: 2876 $

assert (max(cartCells) <= prod(G.cartDims));
assert (min(cartCells) > 0);

actnum                   = false([prod(G.cartDims), 1]);
actnum(G.cells.indexMap) = true;

c           = cumsum(actnum) .* actnum;
activeCells = c(cartCells);
activeCells = activeCells(activeCells > 0);
   
%map = @(i) sparse(1:sum(a(i)), i(a(i)),1, sum(a(i)), numel(a))*cumsum(a);
