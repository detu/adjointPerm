function S = rmSource(S, cells) %#ok
% rmSource -- Remove specified sources/sinks in given system
%
% SYNOPSIS:
%   S = rmSource(S, cells)
%
% PARAMETERS:
%   S     - System structure.
%
%   cells - Vector of indices to cells.
%
% RETURNS:
%   S - System structure with modified S.RHS.g_src, i.e.,
%       S.RHS.g_src(cells) = zeros(size(cells)).
%
% SEE ALSO:
%   putSource, putBC, rmBC.

% $Id: rmSource.m 1432 2009-02-23 15:32:56Z ilig $

error('Function ''rmSource'' is obsolete')

%{
if nargout == 0,
   error(id('NoOutput'), ...
         'Use ''S'' structure both for input and output.')
end

numCells = numel(S.RHS.g_src);
if max(cells) <= numCells,
   S.RHS.g_src(cells) = 0.0;
else
   error(id('Input:Inconsistent'), ...
         'Input is inconsistent.');
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['rmSource:', s];
%}
