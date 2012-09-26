function S = rmBC(S, faces, type) %#ok
% rmBC -- Remove flux or pressure boundary conditions
%
% SYNOPSIS:
%   S = putBC(S, faces, type);
%
% PARAMETERS:
%   S     - System structure.
%
%   faces - Vector of indices to exterior faces for which to remove
%           previously imposed boundary conditions.
%
%   type  - 'pressure' or 'flux'.
%
% RETURNS:
%   S - System structure with modified S.RHS.f_bc or S.RHS.h_bc
%       depending on input 'type'.
%
% SEE ALSO:
%   putBC, putSource, rmSource.

% $Id: rmBC.m 1432 2009-02-23 15:32:56Z ilig $

error('Function ''rmBC'' is obsolete')

%{
if nargout == 0,
   error(id('NoOutput'), ...
         'Use ''S'' structure both for input and output.')
end

S = putBC(S, faces, 0.0, type);
if strcmp(type, 'pressure'),
   S.RHS.dirichletFaces(faces) = false;
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['rmBC:', s];
%}
