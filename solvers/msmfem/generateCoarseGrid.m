function CG = generateCoarseGrid(G, partition, varargin)
%Build coarse grid data structure from partion of existing fine grid.
%
% SYNOPSIS:
%   CG = generateCoarseGrid(G, p)
%   CG = generateCoarseGrid(G, p, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - grid_structure data structure describing fine-scale
%             discretisation of reservoir geometry.
%
%   p       - Partition vector of size [G.cells.num, 1] describing the
%             coarse grid.  We assume that all coarse blocks are connected.
%             The partition vector is typically created by function
%             partitionCartGrid or function partitionUI.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Verbose -- Whether or not to time the coarse grid
%                            generation.
%                            Logical. Default value = FALSE.
%
% RETURNS:
%   CG - Coarse grid structure having the following fields:
%         - cells
%             - num       : Number of coarse blocks.
%             - subCells  : Sparse matrix of size
%                           G.cells.num-by-CG.cells.num defined such that
%                           subCells(i,j)==1 if fine-cell 'i' is part of
%                           coarse block 'j' (and zero otherwise).
%
%         - faces
%             - num       : Number of coarse faces.
%             - neighbors : Coarse blocks sharing common coarse faces.
%                           Size CG.faces.num-by-2.  Specifically, the
%                           coarse blocks 'neighbors(i,1)' and
%                           'neighbors(i,2)' share coarse face 'i'.
%             - tag       : Coarse face tags (inherited from fine grid).
%
%         - cellFaces     : An m-by-2 array of (block)-to-(coarse face)
%                           connections.  Specifically, if
%                           cellFaces(i,1)==j and cellFaces(i,2)==k, then
%                           coarse block 'j' is connected to coarse face
%                           'k'.
%
% SEE ALSO:
%   generateCoarseSystem, generateCoarseWellSystem.

%{
#COPYRIGHT#
%}

% $Id: generateCoarseGrid.m 2933 2009-10-01 15:38:47Z bska $

   opt = struct('Verbose', false);
   opt = merge_options(opt, varargin{:});

   if opt.Verbose, fprintf('Generating coarse grid ... '), tic, end

   %-----------------------------------------------------------------------
   %- Re-index coarse blocks according to pos volume ----------------------
   %
   [partition, numCC] = get_active(partition);

   %-----------------------------------------------------------------------
   %- Detect neighbouring coarse blocks -----------------------------------
   %
   [i, j] = connection(G, partition, numCC);

   %-----------------------------------------------------------------------
   %- Find coarse boundary faces ------------------------------------------
   %
   [ctags, ibdry, jbdry] = boundary_faces(G, partition);

   % Tag all internal coarse faces with tag zero.
   ctags = vertcat(zeros([numel(j), 1]), ctags{:});

   % Add coarse boundary faces to neighbour matrix.
   i = [i; ibdry];
   j = [j; jbdry];

   %-----------------------------------------------------------------------
   %- Build final coarse grid structure -----------------------------------
   %
   CG.cells.num       = numCC;
   CG.cells.subCells  = sub_cells(G, partition, numCC);

   CG.faces.num       = numel(i);
   CG.faces.neighbors = uint32([i, j]);
   CG.faces.tag       = ctags;

   CG.cellFaces       = int32(cell_faces(i, j));
   CG.partition       = partition;

   tocif(opt.Verbose)
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [p, nb] = get_active(p)
   nb = max(p);

   B    = find(accumarray(p, 1));
   I    = zeros([nb, 1]);
   I(B) = 1 : numel(B);

   p  = I(p);
   nb = max(p);      % Number of coarse blocks.
end

%--------------------------------------------------------------------------

function [i, j] = connection(g, p, nb)
   p = [1; p(:) + 1];
   N = p(g.faces.neighbors + 1);
   N = N(N(:,1) ~= N(:,2), :);  % Exclude block-internal connections.

   % Build (symmetric) connectivity matrix.
   conn = sparse([N(:,1); N(:,2); (1 : nb).'], ...
                 [N(:,2); N(:,1); (1 : nb).'], 1);

   % Extract coarse block neighbours to be entered into first ('i') and
   % second ('j') column of array 'CG.faces.neighbors'.
   [j, i] = find(tril(conn(2:end,2:end), -1));
end

%--------------------------------------------------------------------------

function [ctags, i, j] = boundary_faces(g, p)
   % Extract (fine-scale) tags from all exterior (fine-scale) faces and the
   % coarse blocks to which these fine-scale faces are connected.
   %
   ext        = any(g.faces.neighbors == 0, 2);
   block      = zeros([g.faces.num, 1]);
   block(ext) = p(sum(g.faces.neighbors(ext,:), 2));

   f = g.cellFaces(ext(g.cellFaces(:,1)), :);

   % Determine tags on exterior coarse faces.
   %
   %  1) Extract unique (fine-scale) tags for all coarse blocks in 'block'.
   %     Result 'ctags' is a cell array of tags (see ACCUMARRAY).
   %
   ctags = accumarray(block(f(:,1)), f(:,2), [], @(t) { unique(t) });

   %  2) Keep only the blocks for which any tags are assigned
   %     (i.e., the boundary blocks (stored in 'keep')).
   %
   keep  = cellfun(@(t) ~isempty(t), ctags);
   ctags = ctags(keep);

   % Count number of additional basis functions (i.e., faces) on boundary
   % blocks.
   %
   numBF = cellfun(@(t) numel(t), ctags);

   % Create cell-neighbour entries for boundary faces.  The blocks in
   % 'keep' connect to the outside (block zero).  Array 'i' will be entered
   % into the first column of 'CG.faces.neighbors' while array 'j' will be
   % entered into the second column.
   %
   i = zeros   ([sum(numBF),     1]);
   j = rldecode( find(keep), numBF );
end

%--------------------------------------------------------------------------

function subC = sub_cells(g, p, nb)
   subC = logical(sparse(1 : g.cells.num, double(p), 1, g.cells.num, nb));
end

%--------------------------------------------------------------------------

function cf = cell_faces(i, j)
   nn = numel(i);

   %-----------------------------------------------------------------------
   % List local faces and corresponding global index ----------------------
   %
   V = sortrows([ [i, (1 : nn).'];
                  [j, (1 : nn).'] ]);

   i1 = find(V(:,1), 1, 'first'); % First non-zero cell (block).

   cf = V(i1:end,:);
end
