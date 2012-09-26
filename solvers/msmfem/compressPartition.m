function p = compressPartition(p)
%Renumber coarse block partitioning to remove any empty coarse blocks.
%
% SYNOPSIS:
%   p = compressPartition(p)
%
% PARAMETERS:
%   p - Original partition vector, may contain empty coarse blocks.
%
% RETURNS:
%   p - Updated partition vector.
%       Renumbered so as to remove empty coarse blocks.
%
% NOTE:
%   If the original partition does not contain any empty coarse blocks,
%   applying this function is an expensive way of doing nothing.
%
% EXAMPLE:
%   p = partitionCartGrid([4, 4, 1], [2, 2, 2]);
%   [p, compressPartition(p)]
%
% SEE ALSO:
%   partitionCartGrid, partitionUI, processPartition, generateCoarseGrid.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-07 12:48:18 +0200 (ma, 07 sep 2009) $
% $Revision: 2667 $

   if any(p(:) < 1),
      error(msgid('Partition:NonPositive'), ...
           ['Partition vector ''p'' must contain positive integer', ...
            ' entries exclusively']);
   end

   active        = find(accumarray(p(:), 1) > 0);
   compr         = zeros([max(p), 1]);
   compr(active) = 1 : numel(active);

   p = compr(p);
end
