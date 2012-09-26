function S = putBC(S, faces, values, type) %#ok
% putBC -- Put flux or pressure boundary conditions on given system
%
% SYNOPSIS:
%   S = putBC(S, faces, values, type)
%
% PARAMETERS:
%   S      - System structure.
%
%   faces  - Vector of indices to exterior faces.
%
%   values - Vector of  pressure/flux values for faces.
%
%   type   - 'pressure' or 'flux'.
%
% RETURNS:
%   S - System structure with modified S.RHS.f or S.RHS.h
%       corresponding to input.
%
% SEE ALSO:
%   rmBC, putSource, rmSource.

% $Id: putBC.m 1432 2009-02-23 15:32:56Z ilig $

error('function ''putBC'' is obsolete');

%{
if nargout == 0,
   error(id('NoOutput'), ...
         'Use ''S'' structure both for input and output.')
end

types  = { 'pressure', 'flux' };

faces  = reshape(faces , [], 1);
values = reshape(values, [], 1);
typeN  = find(strcmpi(type, types));

if ~isempty(typeN),
   numFaces = numel(S.RHS.h_bc);
   if (max(faces) <= numFaces)         && ...
      ((numel(faces) == numel(values)) || ...
       (numel(values) == 1)),

      [cfIx, trash] = find(S.D(:, faces));    %#ok

      if numel(cfIx) == numel(faces),
         if typeN == 1,
            % Pressure/Dirichlet BC

            S.RHS.f_bc          (cfIx ) = - values;
            S.RHS.dirichletFaces(faces) = true;
            S.RHS.neumannFaces  (faces) = false;
         else
            % Flux/Neumann BC

            S.RHS.h_bc        (faces) = values;
            S.RHS.neumannFaces(faces) = true;
         end
      else
         error(id('Faces:NotBoundary'), ...
               'Given ''faces'' are not on boundary.')
      end
   else
      error(id('Input:Inconsistent'), ...
            'Input is inconsistent.');
   end
else
   error(id('BCType:Unknown'), ...
         'Unknown BC type: ''%s''', type);
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['putBC:', s];
%}
