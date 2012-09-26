function p = partitionCartGrid(cartDims, partDims)
%Partition a Cartesian grid.
%
% SYNOPSIS:
%   partition = partitionCartGrid(cartDims, partDims)
%
% PARAMETERS:
%   cartDims - [nx ny nz] vextor of fine-grid cell dimensions.
%
%   partDims - [cnx cny cnz] vector of coarse-grid block dimensions.
%
% RETURNS:
%   partition - Vector of size [nx*ny*nz 1] with entries equal to
%               coarse block index.
%
% SEE ALSO:
%   processPartition.

%{
#COPYRIGHT#
%}

% $Id: partitionCartGrid.m 1953 2009-03-31 10:54:12Z bska $

nx  = cartDims(1); ny  = cartDims(2); nz  = cartDims(3);
cnx = partDims(1); cny = partDims(2); cnz = partDims(3);

xC = ceil(reshape(repmat(1:nx, 1    , ny*nz), [], 1) ./ (nx / cnx));
yC = ceil(reshape(repmat(1:ny, nx   ,    nz), [], 1) ./ (ny / cny));
zC = ceil(reshape(repmat(1:nz, nx*ny, 1    ), [], 1) ./ (nz / cnz));

p = xC + cnx*((yC-1) + cny*(zC-1));
