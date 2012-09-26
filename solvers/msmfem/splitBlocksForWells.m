function p = splitBlocksForWells(p, splitFact, G, W)
%Split coarse blocks near wells.
%
% SYNOPSIS:
%   p = splitBlocksForWells(p, splitFact, G, W)
%
% PARAMETERS:
%   p         - Original partition vector
%
%   splitFact - Vector of coarse block splitting factors in each physical
%               direction.  Each block determined to be close to a well
%               will be split into PROD(splitFact) new coarse blocks.
%
%   G         - Grid structure as described by grid_structure.
%
%   W         - Well structure defined by addWell &c.
%
% RETURNS:
%   p - Updated partition vector.
%
% NOTE:
%   This function potentially introduces a lot of new coarse blocks which
%   will impede execution speed if subsequently used to define a coarse
%   grid on which the multiscale pressure solver is employed.
%
% SEE ALSO:
%   splitBlocks, refinePartitionForWells.

%{
#COPYRIGHT#
%}

% $Id: splitBlocksForWells.m 1953 2009-03-31 10:54:12Z bska $

p  = [0; p];
cN = p(G.faces.neighbors + 1);
cN = unique(cN(prod(double(cN),2) > 0, :), 'rows');
nb = max(cN(:));
cN = [cN; cN(:,[2,1]); (1 : nb).', (1 : nb).'];
Gn = spones(sparse(cN(:,1), cN(:,2), 1));

wb    = unique(p(vertcat(W.cells) + 1));
v     = zeros([nb, 1]);
v(wb) = 1;

% Neighbourhood == BFS(2) from wb
v = Gn * v; v = Gn * v;

% Only split blocks connected to at least two direct neighbours of 'wb'.
p = splitBlocks(find(v >= 2), p(2:end), splitFact, G);
