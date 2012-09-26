function S = putSource(S, varargin) %#ok
% putSource -- Put source/sink in given system
%
% SYNOPSIS:
%   S = putSource(S, src)
%   S = putSource(S, cells, values)
%
% PARAMETERS:
%   S      - System structure as defined by (e.g.) assembleMimeticSystem.
%
%   And then either
%      src    - Source structure as defined by function 'addSource'.
%
%   Or
%      cells  - Vector of indices to cells.
%
%      values - Vector of source values for given cells.
%
% RETURNS:
%   S - System structure with modified S.RHS.g_src, i.e.,
%       S.RHS.g_src(cells) = values.
%
% SEE ALSO:
%   addSource, rmSource, putBC, rmBC.

% $Id: putSource.m 1432 2009-02-23 15:32:56Z ilig $

error('Function ''putSource'' is obsolete')

%{
if nargout == 0,
   error(id('NoOutput'), ...
         'Use ''S'' structure both for input and output.')
end

if nargin == 2 &&           ...
   isstruct(varargin{1}) && ...
   isfield(varargin{1}, 'cell'),

   src = varargin{1};

   s = accumarray(src.cell, src.rate);
   i = accumarray(src.cell, 1);

   S.RHS.g_src(i > 0) = s(i > 0);
else
   cells  = varargin{1};
   values = varargin{2};

   numCells = numel(S.RHS.g_src);
   if (max(cells) <= numCells) && ...
      ((numel(cells) == numel(values)) || (numel(values) == 1)),
      S.RHS.g_src(cells) = values;
   else
      error(id('Input:Inconsistent'), ...
           'Input is inconsistent.');
   end
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['putSource:', s];
%}
