function G = buildCoordGrid(grdecl)
% buildCoordGrid -- Construct SAMSIM grid from *matching* GRDECL
%                      grid specification.
%
% SYNOPSIS:
%   G = buildCoordGrid(grdecl)
%
% PARAMETERS:
%   grdecl  - Eclipse file output structure as defined by readGRDECL.
%             Must contain at least the fields 'cartDims', 'COORD[XYZ]'
%
% RETURNS:
%   G       - Grid structure as described in grid_structure, though
%             without all of the required geometric information (e.g. cell
%             volumes and cell centroids).
%
% EXAMPLE:
%   casefile = [DATADIR, filesep, 'case.grdecl']
%   grdecl   = readGRDECL(casefile)
%   G        = buildCoordGrid(grdecl)
%   G        = computeGeometry(G)
%
% SEE ALSO:
%   readGRDECL, computeGeometry, grid_structure.

% $Id: buildCoordGrid.m 1780 2009-03-19 13:06:09Z bska $

if ~isfield(grdecl, 'ACTNUM')
   ACTNUM = ones([prod(grdecl.cartDims), 1]);
else
   ACTNUM = grdecl.ACTNUM;
end

%%% Find cells and faces
[C, ncell, N]        = cellCoordinates(grdecl);
[globalFaces, activeFaces, tag] = buildFaces(N);

%%% Build mapping from faces to cells and vice versa
protoCellFaces = kron((1 : ncell) .', ones([6, 1]));
numCellFaces   = accumarray(protoCellFaces(activeFaces), ...
                            1, [ncell, 1]);

%%% Find active cells and faces
active         = ACTNUM & (numCellFaces > 2);
numActiveCells = sum(active);
activeCF       = logical(kron(active, ones([6, 1])) .* activeFaces);
activeCellNo         = zeros([ncell, 1]);
activeCellNo(active) = (1 : numActiveCells) .';

centroidF = C(globalFaces(activeCF,1),:) + ...
            C(globalFaces(activeCF,2),:) + ...
            C(globalFaces(activeCF,3),:) + ...
            C(globalFaces(activeCF,4),:);

%%% Make the face structure
[c, i1, i2] = unique(single(centroidF), 'rows');
ix = find(activeCF);
cellFaces = [activeCellNo(protoCellFaces(activeCF)), ...
             i2,                                     ...
             1-2*mod(ix,2)];

globalFaceNodes = globalFaces(ix(i1), :);
tag = tag(ix(i1));
numFaces  = numel(i1);
faceNodes = [kron((1 : numFaces).', ones([4, 1])), ...
             reshape(globalFaceNodes.', [], 1),    ...
             reshape(globalFaceNodes(:, [2, 3, 4, 1]).', [], 1)];

%%% Assign data to G structure
G.nodes = struct('num', size(C,1), ...
                 'coords', C);
[IJK{1:3}] = ind2sub(grdecl.cartDims(:) .', (1 : ncell) .'); IJK = [IJK{:}];
G.cells = struct('num',      numActiveCells,       ...
                 'numFaces', uint32(numCellFaces(active)), ...
                 'indexMap', uint32(find(active)),         ...
                 'ijkMap',   IJK(active,:));
I1 = find(cellFaces(:,3) > 0);
I2 = find(cellFaces(:,3) < 0);
neighbors                     = zeros([numFaces, 2]);
neighbors(cellFaces(I1,2), 1) = cellFaces(I1,1);
neighbors(cellFaces(I2,2), 2) = cellFaces(I2,1);

G.faces = struct('num',       numFaces,                 ...
                 'numNodes',  repmat(uint8(4), [numFaces, 1]), ...
                 'neighbors', uint32(neighbors),...
                 'tag',       uint8(tag));

G.cellFaces =  int32(cellFaces(:,2:3));
G.faceNodes = uint32(faceNodes(:,2:3));
G.cartDims  = grdecl.cartDims(:) .';
G.type = 'buildCoordGrid';

%--------------------------------------------------------------------------
% Private helpers follow
%--------------------------------------------------------------------------


function [C, ncell, nodeIndex] = cellCoordinates(grdecl)
nx = grdecl.cartDims(1);
ny = grdecl.cartDims(2);
nz = grdecl.cartDims(3);

C = [ grdecl.COORDX(:), grdecl.COORDY(:), grdecl.COORDZ(:)];

[C,I,J] = unique(C(:,[3 2 1]),'rows'); %#ok
C = C(:,[3 2 1]);
nodeNo = reshape(J,grdecl.cartDims+1);

n1 = nodeNo(1:nx,   1:ny,   1:nz  );
n2 = nodeNo(2:nx+1, 1:ny,   1:nz  );
n3 = nodeNo(1:nx,   2:ny+1, 1:nz  );
n4 = nodeNo(2:nx+1, 2:ny+1, 1:nz  );
n5 = nodeNo(1:nx,   1:ny,   2:nz+1);
n6 = nodeNo(2:nx+1, 1:ny,   2:nz+1);
n7 = nodeNo(1:nx,   2:ny+1, 2:nz+1);
n8 = nodeNo(2:nx+1, 2:ny+1, 2:nz+1);

nodeIndex = [n1(:), n2(:), n3(:), n4(:), ...
             n5(:), n6(:), n7(:), n8(:)];

ncell = prod(grdecl.cartDims);

%--------------------------------------------------------------------------
function [F, activeF, tag] = buildFaces(N)
F = reshape(N(:,[1, 3, 7, 5,   2, 4, 8, 6,...  
                 1, 5, 6, 2,   3, 7, 8, 4,... 
                 1, 2, 4, 3,   5, 6, 8, 7])', 4, [])';
tag = repmat([1:6]',size(N,1),1); %#ok
activeF = cellfun( @(a)( numel(unique(a)) ),...
                   mat2cell(F,ones(size(F,1),1),4) ) > 2;

