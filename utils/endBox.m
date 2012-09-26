function c = endBox(cartDims)
%Restore input box to default settings for given Cartesian size.
%
% SYNOPSIS:
%   cells = endBox(cartDims)
%
% PARAMETERS:
%   cartDims - Cartesian dimensions of global index space.  Assumed to be a
%              d-element vector of positive integers.  This parameter often
%              corresponds to the [nx, ny, nz] triplet of the ECLIPSE
%              'SPECGRID' (or 'DIMENS') keyword.
%
% RETURNS:
%   cells - List of linear indices corresponding to entire Cartesian box.
%
% SEE ALSO:
%   gridBox.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-29 13:22:23 +0200 (ti, 29 sep 2009) $
% $Revision: 2889 $

   c = (1 : prod(double(cartDims))) .';
end
