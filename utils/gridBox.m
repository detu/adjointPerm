function c = gridBox(cartDims, region)
%Construct list of cells corresponding to index space subset.
%
% SYNOPSIS:
%   cells = gridBox(cartDims, region)
%
% PARAMETERS:
%   cartDims - Cartesian dimensions of global index space.  Assumed to be a
%              d-element vector of positive integers.  This parameter often
%              corresponds to the [nx, ny, nz] triplet of the ECLIPSE
%              'SPECGRID' (or 'DIMENS') keyword.
%
%   region   - Array of lower/upper bounds on Cartesian index space subset
%              of 'cartDims'.  Assumed to be a NUMEL == 2*d array of
%              positive integers whose consecutive linear indices
%              correspond to separate index space dimensions.  In
%              particular, the indices
%
%                      region(2*(i-1) + 1) : region(2*i)
%
%              define the subset of the i'th Cartesian dimension.
%
% RETURNS:
%   cells - List of linear indices corresponding to the 'region' subset.
%           This return value is suitable for passing as the second
%           parameter to function 'cart2active' in order to extract a
%           subset of *active* cells from a corner-point grid data
%           structure as defined by 'grid_structure'.
%
% EXAMPLE:
%   % Extract inner 3-by-3-by-3 cells of 5-by-5-by-5 block.
%   region = repmat([2, 4], [1, 3])     % Lower = 2, Upper = 4.
%   cells  = gridBox([5, 5, 5], region)
%
% NOTE:
%   Function gridBox supports 'cartDims' of general dimension.
%
% SEE ALSO:
%   endBox, cart2active.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-29 13:22:23 +0200 (ti, 29 sep 2009) $
% $Revision: 2889 $

   d = numel(cartDims);

   % Assert correct dimensions.
   assert (d > 0);
   assert (numel(region) == 2*d);

   cartDims = reshape(cartDims, [], 1);
   region   = reshape(region  , 2, []) .';

   % Assert validity of subset.
   assert (all(all(1 <= region)));
   assert (all(all(bsxfun(@le, region, cartDims))));

   % Assert increasing index ranges.
   assert (all(diff(region, [], 2) >= 0));

   % Compute valid subset.
   x = arrayfun(@colon, region(:,1), region(:,2), ...
                'UniformOutput', false);
   [ix{1:d}] = ndgrid(x{:});
   c = reshape(sub2ind(cartDims .', ix{:}), [], 1);
end
