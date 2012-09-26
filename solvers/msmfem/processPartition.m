function partition2 = processPartition(G, partition, varargin)
%Split disconnected coarse blocks into new blocks.
%
% SYNOPSIS:
%   p2 = processPartition(G, p)
%   p2 = processPartition(G, p, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   p       - Vector of size G.cells.num-by-1 of initial cell-to-block
%             mappings.  Assumed to be a vector of positive integers.  If
%             there are no non-connected coarse blocks, then ALL(p2 == p),
%             meaning the coarse blocks numbers are preserved for all
%             internally connected coarse blocks.
%
%   'pn'/pv - List of 'key'/value pairs designating optional parameters.
%             Currently supported parameters are
%               - Verbose   -- Whether or not to emit progress reports
%                              during the computation.
%                              Logical.  Default value: FALSE.
%
%               - Reconnect -- Whether or not to perform a 'reconnection'
%                              step in which cells that have not been
%                              assigned a new (coarse) block will be
%                              connected to the nearest block.
%                              Logical.  Default value: TRUE.
%
%                              If Reconnect==FALSE, (i.e., don't reconnect
%                              disconnected cells) each disconnected cell
%                              will be assigned a separate new block.
%
% RETURNS:
%   p2 - Updated partition with only connected coarse blocks.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-20 18:22:24 +0200 (ti, 20 okt 2009) $
% $Revision: 3042 $

opt = struct('Verbose', false, 'Reconnect', true);
opt = merge_options(opt, varargin{:});

if any(accumarray(partition, 1) == 0),
   error(msgid('Partition:EmptyBlock'), ...
        ['At least one of the coarse blocks does not ', ...
         'contain any cells.  Use ''compressPartition''.']);
end

partition = partition - min(partition(:)) + 1;

%% Establish coarse block internal topology
[N, blockFaces, bfIx] = getBlockNeighbours(G, partition);

%% Main processing loop
% Objective: Split non-connected coarse blocks into new coarse blocks while
% preservering the original coarse block number for connected coarse
% blocks.
%
nBlk   = max(partition);
maxBlk = nBlk;

if opt.Verbose, h = waitbar(0, 'Processing partition...'); tic, end

partition2 = repmat(-1, size(partition));
for b = 1 : nBlk,
   %  1) Extract all internal faces in block 'b'.
   %
   bf = blockFaces(bfIx(b) + 1 : bfIx(b + 1));

   %  2) Extract the specific cells currently in block 'b'.
   %
   bN  = N(bf,:);                    % Neighbouring cells within block.
   bc  = unique(reshape(bN, [], 1)); % Only unique cells for renumbering.
   nbc = numel(bc);

   %  3) Map global cell numbers to block-local cell numbers and
   %     re-connect the global faces within the block (i.e., 'bf') to the
   %     local cells to maintain block-internal connectivity (i.e.,
   %     topology).
   %
   g2l = sparse(double(bc), 1, 1 : nbc);
   bf  = reshape(full(g2l(bN)), size(bN));

   %  4) Construct symmetric adjacency matrix for block 'b' (non-zero
   %     diagonal) connectivity.
   % Ref:
   %    <URL:http://blogs.mathworks.com/steve/2007/03/20/
   %         connected-component-labeling-part-3/#comments>,
   %    specifically comment No. 8.
   %
   adj = sparse([bf(:,1); bf(:,2); (1 : nbc)'], ...
                [bf(:,2); bf(:,1); (1 : nbc)'],  1 );

   %  5) Compute the connected components in this coarse block by means
   %     of Dulmage-Mendelsohn permutation.
   %
   [p, r, r] = dmperm(adj);   nComp = numel(r) - 1;   % [p,r,r] is no typo!

   %  6) Assign new block numbers to the individual connected components
   %     found in step 5 while taking care to preserve the original block
   %     number of the first (and, possibly, only) block (i.e., component).
   %     This ensures that the output partition equals the input partition
   %     if no coarse blocks must be split.  As a special case, array
   %     'blkNo' must be empty when nbc==0 to avoid triggering an assertion
   %     in 'rldecode'.  This condition occurs if coarse block 'b' contains
   %     a single cell.
   %
   blkNo             = [b(nbc > 0), maxBlk + (1 : nComp - 1)];
   partition2(bc(p)) = rldecode(blkNo, reshape(diff(r), 1, []), 2);

   maxBlk = maxBlk + max(nComp - 1, 0);

   if opt.Verbose, waitbar(b / nBlk, h), end
end

if opt.Verbose, toc, close(h); end

%--------------------------------------------------------------------------
% Handle single-cell coarse blocks in original partition ------------------
%
% Single cells already constituting separate coarse blocks should remain so
% when processed.  However, the above algorithm does not touch such cells
% as there are no *internal* faces within single-cell coarse blocks.  Take
% care to preserve the original coarse block numbers of these single cell
% coarse blocks.
%
single_p                                = false([nBlk, 1]);
single_p(accumarray(partition, 1) == 1) = true;
single_p                                = single_p(partition);
partition2(single_p)                    = partition(single_p);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Reconnect any cells yet to be repartitioned (due to e.g. lacking contacts
% to other cells in their respective original coarse blocks) to appropriate
% new coarse blocks.
%
noblk = find(partition2 == -1);

if opt.Reconnect,
   k = 0;
   while ~isempty(noblk),
      if k > 0,
         dispif(opt.Verbose, 'numel(noblk) == %d\n', numel(noblk));
      end
      partition2 = reconnect(G, partition2, noblk);
      noblk = find(partition2 == -1);
      k = k + 1;
   end
else
   partition2(noblk) = max(partition2) + (1 : numel(noblk));
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function [N, blockFaces, bfIx] = getBlockNeighbours(G, partition)
N = G.faces.neighbors(all(G.faces.neighbors > 0, 2), :);

% Extract faces entirely internal to a block.
nblk = max(partition);
pN = partition(N);
iF = pN(:,1) == pN(:,2);
N  = N (iF,:);
pN = pN(iF,:);

% Split faces (contacts) according to coarse block numbers.
% The block faces in [bfIx(b) + 1 : bfIx(b+1)] are all internal faces in
% block 'b'.
%
blockFaces = sortrows([pN(:,1), (1 : size(pN,1)) .']);
bfIx       = cumsum([0; accumarray(blockFaces(:,1), 1, [nblk, 1])]);
blockFaces = blockFaces(:,2);

%--------------------------------------------------------------------------

function partition2 = reconnect(G, partition2, noblk)
noblk = reshape(noblk, 1, []);  % Row vector for direct FOR loop indexing
for c = noblk,
   % Determine all neighbouring cells of cell 'c' (-> 'other').
   %
   faces = G.cellFaces(G.cells.facePos(c) : G.cellsfacePos(c + 1) - 1, 1);
   other = reshape(G.faces.neighbors(faces,:), [], 1);

   % Exclude 'other' cells which are either external (0), equal to myself
   % (c) or not a member of any coarse blocks themselves (partition2==-1).
   %
   other(ismember(other, [0; c])) = [];
   other(partition2(other) == -1) = [];

   % Connect cell 'c' to first 'other' cell (if any other cells remain).
   %
   if ~isempty(other), partition2(c) = partition2(other(1)); end
end
