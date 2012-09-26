function bf = extractWellBF(basis, sz)
%Form matrix of well basis function values.
%
% SYNOPSIS:
%   bf = extractWellBF(basis, m)
%
% PARAMETERS:
%   basis   - Cell array of packed basis function values as defined by
%             function 'evalWellBasis'.  May be either of the flux or
%             pressure basis functions (well fields 'CS.basis' or
%             'CS.basisP').  The cell array is assumed to be restricted to
%             only those basis functions which will be entered into the
%             resulting matrix.
%
%   m       - Number of rows in the resulting (sparse) basis function value
%             matrix.  This number differs if we're constructing flux
%             values (m==number of fine-scale half-contacts) or pressure
%             values (m==number of fine-scale cells).
%
% RETURNS:
%   bf - Matrix, size m-by-NUMEL(basis), of respective basis function
%        values.  The values of the j'th basis function (i.e., basis{j}) is
%        stored in bf(:,j).
%
% NOTE:
%   This function only expands the packed storage format represented by the
%   cell array 'basis'.  This impacts the return values of 'extractBF' when
%   used to derived flux basis function values from the 'CS.basis' field.
%   As function 'evalWellBasis' returns flux basis function values in the
%   form 'B*v', these are the values which are entered into the 'bf'
%   matrix.  It is the caller's responsibility to compute the required
%   matrix product, e.g., S.BI * bf, in order to derive the correct flux
%   basis function values.
%
% SEE ALSO:
%   extractBF.

%{
#COPYRIGHT#
%}

% $Id: extractWellBF.m 1953 2009-03-31 10:54:12Z bska $

   [i,j,v] = cellfun(@(c) c{1:3}, basis, 'UniformOutput', false);
   bf      = sparse(       vertcat(i{:}) , ...
                    cumsum(vertcat(j{:})), ...
                           vertcat(v{:}) , sz, numel(j));
end
