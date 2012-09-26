function blockIx = partitionUI(G, coarseDim)
%Partition grid uniformly in logical space.
%
% SYNOPSIS:
%   blockIx = partitionUI(G, coarseDim)
%
% DESCRIPTION:
%   Partitions a corner point grid in grid_structure format (relatively)
%   uniformly in index (logical) space.  Each coarse block is
%   comprised of roughly the same number of cells from the fine grid.
%
% PARAMETERS:
%   G         - grid_structure structure having either a valid field
%               'G.cells.ijkMap' or both valid fields 'G.cartDims' and
%               'G.cells.indexMap'.
%
%   coarseDim - Number of coarse blocks in each physical direction.
%               Assumed to be a LENGTH 2 or 3 vector of integers.
%
% RETURNS:
%   blockIx   - A G.cells.num-by-1 vector mapping cells to coarse blocks.
%
% SEE ALSO:
%   processPartition, partitionLayers.

%{
#COPYRIGHT#
%}

% $Date: 2009-08-27 09:06:36 +0200 (to, 27 aug 2009) $
% $Revision: 2599 $

if ~grid_ok(G),
   error(msgid('grid_structure:Invalid'), ...
         'Grid is not a valid grid_structure structure');
end

if ~any(numel(coarseDim) == [2, 3]) || any(coarseDim < 1),
   error(msgid('coarseDim:Invalid'), ...
        ['Parameter ''coarseDim'' must be a two or three component ', ...
         'vector\ncontaining strictly positive integers.']);
end

if isfield(G.cells, 'ijkMap');
   ijk = double(G.cells.ijkMap);
   M   = double(max(ijk));
else
   [ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
   ijk = double([ijk{:}]);
   M = ([max(ijk(:,1))-min(ijk(:,1))+1, max(ijk(:,2))-min(ijk(:,2))+1, ...
         max(ijk(:,3))-min(ijk(:,3))+1]);
end

blockIx = zeros([G.cells.num, 1]);
for d = numel(coarseDim) : -1 : 1,
   B = coarseDim(d);
   blockIx = lbLinDist(ijk(:,d) - min(ijk(:,d)), M(d), B) + B*blockIx;
end
blockIx = blockIx + 1;  % Map block numbers to (1..PROD(coarseDim))


function f = lbLinDist(f, M, B)
% lbLinDist -- Load-balanced linear distribution
% See Eric F. Van de Velde, Concurrent Scientific Computing,
% 1994, Springer Verlag, p. 54 (Sect. 2.3) for details.
%
% Maps index set (0..M-1) to blocks (0..B-1).

L = floor(M ./ B);  % Tentative number of cells per coarse block.
R = mod(M, B);      % Additional cells not previously accounted for.
f = max(floor(f ./ (L + 1)), floor((f - R) ./ L));


function bool = grid_ok(G)
bool = ~isempty(G)          && isstruct(G)       && ...
        isfield(G, 'cells') &&                      ...
       ~isempty(G.cells)    && isstruct(G.cells);
bool = bool && ((isfield(G.cells, 'ijkMap') && ...
                 size(G.cells.ijkMap,1) == G.cells.num) || ...
                (isfield(G, 'cartDims') && isfield(G.cells, 'indexMap')));
