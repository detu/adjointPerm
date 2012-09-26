function CS = assignBasisFuncs(CS, V, P)
%Update reservoir coarse system basis functions following 'evalBasisFunc'.
%
% SYNOPSIS:
%   CS = assignBasisFuncs(CS, V, P)
%
% PARAMETERS:
%   CS   - Original coarse system.  Must have allocated fields '.basis'
%          and '.basisP' of correct sizes.
%
%   V, P - New basis function values for flux (V) and pressure (P).
%          Assumed to be the return values from function 'evalBasisFunc'.
%
% RETURNS:
%   CS - Updated reservoir coarse system structure with new values for flux
%        and pressure basis functions.
%
% SEE ALSO:
%   generateCoarseSystem, evalBasisFunc.

%{
#COPYRIGHT#
%}

% $Id: assignBasisFuncs.m 1953 2009-03-31 10:54:12Z bska $

   % Extract coarse faces for which new basis function values are defined.
   f = cellfun(@(x) x{3}, V);

   CS.basis (f) = V;
   CS.basisP(f) = P;
end
