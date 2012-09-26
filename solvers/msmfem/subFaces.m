function [nsub, sub] = subFaces(g, cg)
%Extract fine-grid faces constituting individual coarse grid faces.
%
% SYNOPSIS:
%   [nsub, sub] = subFaces(G, CG)
%
% PARAMETERS:
%   G  - Grid data structure as described by grid_structure.
%   CG - Coarse grid data structure.
%
% RETURNS:
%   nsub - Number of fine-grid faces belonging to each individual coarse
%          grid face.  Specifically, nsub(i) is the number of fine-grid
%          faces belonging to coarse face 'i'.
%
%   sub  - Actual fine-grid subfaces represented in a condensed storage
%          format.  Specifically, if IX = CUMSUM([0; nsub]), then the
%          subfaces of coarse face 'i' are sub(IX(i) + 1 : IX(i + 1)).
%
% EXAMPLE:
%   % Generate moderately large grid and corresponding coarse grid
%   G  = cartGrid([60, 220, 85]);
%   p  = partitionUI(G, [5, 11, 17]);
%   CG = generateCoarseGrid(G, p);
%
%   % Extract subfaces
%   [nsub, sub] = subFaces(G, CG);
%
%   % Find fine-faces belonging to coarse face 10
%   sub_ix = cumsum([0; nsub]);
%   sf10   = sub(sub_ix(10) + 1 : sub_ix(10 + 1))
%
%   % Find index of coarse face to which fine-face 7 belongs
%   f2c = sparse(sub, 1, rldecode((1 : CG.faces.num) .', nsub));
%   cf7 = f2c(7)
%
% SEE ALSO:
%   grid_structure, generateCoarseGrid.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

%--------------------------------------------------------------------------
% Extract extended partition vector.
%
% Note: We map the outside (cell zero) to coarse block zero too.
%
[i,j]  = find(cg.cells.subCells);
p      = zeros([size(cg.cells.subCells,1) + 1, 1]);
p(i+1) = j;

%--------------------------------------------------------------------------
% Build coarse topology and face tag table.
%
% The table is defined by
%
%   A = [p(N), tag, face]
%
% whence row 'i' contains the coarse block (b1,b2) sharing half-face 'i',
% the tag (t) on (half-face) 'i' and the corresponding global fine-scale
% face of half-face 'i'.  Note specifically that all *internal* faces, i.e.
% those faces shared by two non-zero cells, get face-tag zero.  This is
% because the existing coarse grid (cg) tags all of its internal (coarse)
% faces with tag zero.
%
% Moreover, the global fine-scale face (final column of 'A') is only really
% needed when we form the output 'sub' as 'face' represents the information
% we eventually require.  However, we need to reduce this data in the same
% manner as the (b1,b2,t) triplets so we carry the information along with
% those triplets.
%
A = [p(g.faces.neighbors(g.cellFaces(:,1),:) + 1), g.cellFaces(:,[2,1])];
A(all(A(:,[1,2]) > 0, 2), 3) = 0;

%--------------------------------------------------------------------------
% Extract only unique coarse faces in the underlying fine grid.
%
% The unique coarse faces are defined as:
%   1) Faces shared by exactly two different coarse blocks.  The outside
%      (i.e., coarse block zero) counts as a different coarse block in this
%      respect so we also take care of coarse boundary faces.
%
%   2) Having a unique face tag (only applicable to boundary faces) meaning
%      we can distinguish different parts of an external coarse face.
%
A          = A(A(:,1) ~= A(:,2), :);            % Only faces where b1 ~= b2
A(:,[1,2]) = sort(A(:,[1,2]), 2);               % Deterministic block sort
[Ar, i, i] = unique(double(A(:,1:3)), 'rows');  % Unique (b1,b2,t) triplets

% Check if we've found the correct number of unique faces.
assert (size(Ar,1) == cg.faces.num);

%--------------------------------------------------------------------------
% Locate and enumerate the coarse faces [(b1,b2,t) triplets] of the cg.
%
% Assert that all of these faces are amongst those faces derived from the
% partition and the underlying fine grid.
%
B          = [sort(double(cg.faces.neighbors), 2), double(cg.faces.tag)];
[loc, loc] = ismember(B, Ar, 'rows');  assert (all(loc > 0));

% Renumber the faces in 'Ar' to correspond to those of cg.
%
cfno      = zeros([cg.faces.num, 1]);
cfno(loc) = 1 : numel(cfno);

%--------------------------------------------------------------------------
% Form output arrays.
%
% Algorithm:
%   1) Extract the unique (coarse-face, fine-face) pairs, sorted by
%      coarse-face index (and then by fine-scale face index).
sub  = unique(double([cfno(i), A(:,end)]), 'rows');

%   2) Count the number of fine-scale faces in each coarse face.
nsub = accumarray(sub(:,1), 1);

%   3) Return only those (sorted) fine-scale faces.
sub  = sub(:,2);
