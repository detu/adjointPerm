function [S, CS] = putSourceMS(G, S, CG, CS, rock, varargin) %#ok
% putSourceMS -- Put source/sink in given coarse and fine system
%
% SYNOPSIS:
%   [S, CS] = putSourceMS(G, S, CG, CS, rock, src)
%   [S, CS] = putSourceMS(G, S, CG, CS, rock, cells, values, sat)
%
% PARAMETERS:
%   G      - Fine grid
%
%   S      - Fine system structure
%
%   CG     - Coarse grid
%
%   CS     - Coarse system structure
%
%   rock   - Rock data structure with valid fields 'perm' and 'poros'.
%
%   cells  - Vector of indices to cells (in G)
%
%   values - Vector of source values for given cells
%
%   sat    - 
%
% RETURNS:
%   S  - System structure with modified S.RHS.g_src.
%   CS - System structure with modified CS.RHS.g_src, and modified CS.basis
%
% SEE ALSO:
%   rmSourceMS, putBCMS, rmSource, putBC, rmBC

% $Id: putSourceMS.m 1829 2009-03-23 14:01:07Z bska $

error('Function ''putSourceMS'' is obsolete')

%{
if nargout ~= 2,
   error(id('Output:Inconsistent'), ...
         'Use ''[S, CS]'' structures both for input and output.')
end

src = get_src(varargin{:});

S   = putSource(S, src);
CS.RHS.g_src = CG.cells.subCells.' * S.RHS.g_src;

% Determine weighting function with which to update basis functions:
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

% Determine blocks for which the basis functions must be updated
[trash, blocks] = find(CG.cells.subCells(src.cell,:));
blocks = unique(blocks);

CS = generateBasis(G, S, CG, CS, weight, ...
                   'Verbose', true,      ...
                   'Blocks', blocks);

%--------------------------------------------------------------------------
% Private helpers follow
%--------------------------------------------------------------------------

function s = id(s)
s = ['putSourceMS:', s];

%--------------------------------------------------------------------------

function src = get_src(varargin)
if nargin == 1,
   src = varargin{1};
else
   assert (nargin == 3);
   src = addSource([], varargin{:});
end
%}
