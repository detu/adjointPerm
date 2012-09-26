function s = cellNodes(g)
%Extract local-to-global vertex numbering for grid cells.
%
% SYNOPSIS:
%   cn = cellNodes(G)
%
% PARAMETERS:
%   G - Grid data structure geometrically discretising a reservoir model.
%
% RETURNS:
%   cn - An m-by-3 array mapping cell numbers to vertex numbers.
%        Specifically, if cn(i,1)==j and cn(i,3)==k, then global vertex 'k'
%        is one of the corners of cell 'j'.  The local vertex number of
%        global node 'k' within cell 'j' may be computed using the
%        statements
%
%           n_vert = accumarray(cn(:,1), 1, [G.cells.num, 1]);
%           offset = rldecode(cumsum([0; n_vert(1:end-1)]), n_vert);
%           loc_no = (1 : size(cn,1)).' - offset;
%
%        This calculation is only for illustration purposes.  It is assumed
%        that the local vertex number will be implicitly available in most
%        applications (e.g., finite element methods).
%
%        Alternatively, the local vertex number is found in cn(i,2). ;). In
%        a cornerpoint grid, local vertex number of zero indicates that the
%        vertex is not part of the original 8 vertices.
%
% SEE ALSO:
%   grid_structure.

%{
#COPYRIGHT#
%}

% $Id: cellNodes.m 2321 2009-06-04 07:21:32Z afg $

   dim = numel(g.cartDims);

   cellNo = rldecode(1:g.cells.num, double(g.cells.numFaces), 2) .';
   cf     = double(g.cellFaces);
   ne     = double(g.faces.numNodes);

   % expand table of cell number and halfface info
   a   = rldecode([cellNo, cf(:,2)], ne(cf(:,1)), 1);

   % the index i expands g.faceNodes to 'halffaceNodes'
   pos = cumsum([0; ne]);
   i   = mcolon(pos(cf(:,1)) + 1, pos(cf(:,1) + 1)) .';


   % binary magic applied to halfface-tag to discover nodes shared by
   % three faces.  Remember that sparse adds contributions when
   % (cell, node)-pairs are repeated.
   numbers = 2 .^ (0 : max(a(:,2)) - 1) .';
   mat     = sparse(a(:,1), double(g.faceNodes(i,1)), numbers(a(:,2)));

   [cellNumber, nodeNumber, tag] = find(mat);

if dim==2
   % For Cartesian grid, these are halfface tags for each of eight
   % vertices
   k = [1, 3; ...  % west-south  = 2^0+2^2+2^4
        2, 3; ...  % east-south
        1, 4; ...  % west-north
        2, 4];     % east-north

elseif dim==3
   % For Cartesian grid, these are halfface tags for each of eight
   % vertices
   k = [1, 3, 5; ...  % west-south-down  = 2^0+2^2+2^4
        2, 3, 5; ...  % east-south-down
        1, 4, 5; ...  % west-north-down
        2, 4, 5; ...  % east-north-down
        1, 3, 6; ...  % west-south-up
        2, 3, 6; ...  % east-south-up
        1, 4, 6; ...  % west-north-up
        2, 4, 6];     % east-north-up
end
   localnode = sparse(sum(2 .^ (k-1), 2), 1, 1 : size(k,1));

   s = sortrows([cellNumber(:), full(localnode(tag(:))), nodeNumber(:)]);
   assert (all(s(:,3) > 0));
end
