function [S, CS] = rmSourceMS(G, S, CG, CS,rock, cells) %#ok
% PUTSOURCE -- Remove source/sink in given coarse and fine system
%
% SYNOPSIS:
%   [S, CS] = rmSourceMS(G, S, CG, CS, cells, values, rock);
%
% PARAMETERS:
%   G          - Fine grid
%
%   S          - Fine system structure
%
%   CG         - Coarse grid
%
%   CS         - Coarse system structure
%
%   cells      - Vector of indices to cells (in G)
%
%   values     - Vector of source values for given cells
%
%   rock       - Rock data structure with valid fields 'perm' and 'poros'.
%
% RETURNS:
%   S  - System structure with modified CS.RHS.g_src
%   CS - System structure with modified CS.RHS.g, and modified CS.basis
%
%
% SEE ALSO:
%   putSourceMS, putBCMS, rmSource, putBC, rmBC,

% $Id: rmSourceMS.m 1829 2009-03-23 14:01:07Z bska $

error('Function ''rmSourceMS'' is obsolete')

%{
if nargout ~= 2,
   error(id('NoOutput'), ...
         'Use ''[CS, S]'' structures both for input and output.')
end

S=rmSource(S, cells);
CS.RHS.g_src = CG.cells.subCells.' *S.RHS.g_src; 

% find weight to update basis functions:
ix = strmatch(CS.basisWeighting, {'perm', 'poros', 'unit'}, 'exact');
switch ix,
   case 1,
      if size(rock.perm, 2) == 1,
         weight = rock.perm;
      elseif size(rock.perm, 2) == 3,
         weight = sum(rock.perm, 2);
      else
         weight = sum(rock.perm(:, [1, 4, 6]), 2);
      end
   case 2,
      weight = rock.poro;
   case 3,
      weight = ones([G.cells.num, 1]);
end
%find blocks where basis functions must be updated
[trash, blocks] = find(CG.cells.subCells(cells,:));
blocks = unique(blocks);

CS = generateBasis(G, S, CG, CS, weight,'Verbose', true, 'Blocks', blocks);

function s = id(s)
s = ['rmSourceMS:', s];
%}
