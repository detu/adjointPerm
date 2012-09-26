function src = evalBasisSource(G, weighting, rock)
%Evaluate synthetic basis function driving source term.
%
% SYNOPSIS:
%   src = evalBasisSource(G, w, rock)
%
% PARAMETERS:
%   G    - Grid data structure.
%   w    - Source weighting mode.  May be one of:
%            - Name (string) of particular weighting mode.  Must be one of
%                - 'perm' -- Weigh sources according to TRACE(K).
%                - 'poro'/'poros' --
%                            Weigh sources according to porosity.
%                - 'unit' -- Uniform sources.
%
%            - Array of already evaluated synthetic source terms.  This
%              array will be returned unchanged, and may be useful when
%              updating basis functions to account for compressible
%              effects.  The array is assumed to contain one non-negative
%              scalar value for each cell in the model.
%
%   rock - Rock data structure.  May be empty, but in the case of
%            - w == 'perm' -> must contain valid field 'rock.perm'.
%            - w == 'poro'/'poros'
%                          -> must contain valid field 'rock.poro'.
%
% RETURNS:
%   src  - Synthetic basis function source term for use with generator
%          functions 'evalBasisFunc' and 'evalWellBasis'.
%
% SEE ALSO:
%   generateCoarseSystem, generateCoarseWellSystem.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 17:21:14 +0200 (to, 01 okt 2009) $
% $Revision: 2932 $

   if isnumeric(weighting) && numel(weighting) == G.cells.num,
      src = weighting;
   elseif ischar(weighting),
      ix = strmatch(weighting, {'perm', 'poro', 'poros', 'unit'}, 'exact');

      if isempty(ix),
         error(msgid('BasisWeighting:Unknown'), ...
               'Basis function weighting strategy ''%s'' is unknown', ...
               weighting);
      end

      switch ix,
         case 1,
            if isempty(rock) || ~isfield(rock, 'perm'),
               error(msgid('RockData:Empty'), ...
                    'Rock data must contain valid field ''perm''.');
            end

            dim = size(G.nodes.coords, 2);
            K   = permTensor(rock, dim);
            src = K(:, 1 : dim+1 : end) * ones([dim, 1]);  % == TRACE(K)

         case {2, 3},
            if isempty(rock) || ~isfield(rock, 'poro'),
               error(msgid('RockData:Empty'), ...
                     'Rock data must contain valid field ''poro''.');
            end

            src = rock.poro;

         case 4,

            src = ones([G.cells.num, 1]);

      end
   else
      error(msgid('BasisWeighting:NotSupported'), ...
           ['Basis function weighting does not conform to any ', ...
            'supported mode.']);
   end
end
