function G = buildMatchingGrid(grdecl)
% buildMatchingGrid -- Construct SAMSIM grid from *matching* GRDECL
%                      grid specification.
%
% SYNOPSIS:
%   G = buildMatchingGrid(grdecl)
%
% PARAMETERS:
%   grdecl  - Eclipse file output structure as defined by readGRDECL.
%             Must contain at least the fields 'cartDims', 'COORD' and
%             'ZCORN'.
%
% RETURNS:
%   G       - Grid structure as described in grid_structure, though
%             without all of the required geometric information (e.g. cell
%             volumes and cell centroids).
%
% EXAMPLE:
%   casefile = [DATADIR, filesep, 'case.grdecl']
%   grdecl   = readGRDECL(casefile)
%   G        = buildMatchingGrid(grdecl)
%   G        = computeGeometry(G)
%
% SEE ALSO:
%   readGRDECL, computeGeometry, grid_structure.

% $Id: buildMatchingGrid.m 1780 2009-03-19 13:06:09Z bska $

if ~isfield(grdecl, 'ACTNUM')
   ACTNUM = ones([prod(grdecl.cartDims), 1]);
else
   ACTNUM = grdecl.ACTNUM;
end

[C, ncell, nnod]           = cellCoordinates(grdecl);
[globalFaces, activeFaces] = buildFaces(C, ncell, nnod);

protoCellFaces = kron((1 : ncell) .', ones([6, 1]));
numCellFaces   = accumarray(protoCellFaces(activeFaces), ...
                            1, [ncell, 1]);

active         = ACTNUM & (numCellFaces > 2);
numActiveCells = sum(active);
activeCF       = logical(kron(active, ones([6, 1])) .* activeFaces);

activeCellNo         = zeros([ncell, 1]);
activeCellNo(active) = (1 : numActiveCells) .';

centroidF = C(globalFaces(activeCF,1),:) + ...
            C(globalFaces(activeCF,2),:) + ...
            C(globalFaces(activeCF,3),:) + ...
            C(globalFaces(activeCF,4),:);

[c, i1, i2] = unique(single(centroidF), 'rows');

ix = find(activeCF);
cellFaces = [activeCellNo(protoCellFaces(activeCF)), ...
             i2,                                     ...
             1-2*mod(ix,2)];


globalFaceNodes = globalFaces(ix(i1), :);
numFaces  = numel(i1);
faceNodes = [kron((1 : numFaces).', ones([4, 1])), ...
             reshape(globalFaceNodes.', [], 1),    ...
             reshape(globalFaceNodes(:, [2, 3, 4, 1]).', [], 1)];

activeEdges = sqrt(sum((C(faceNodes(:,2),:) - ...
                        C(faceNodes(:,3),:)) .^ 2, 2)) > sqrt(eps);
faceNodes = faceNodes(activeEdges,:);

G.nodes = struct('num', size(C,1), 'coords', C);

[IJK{1:3}] = ind2sub(grdecl.cartDims(:) .', (1 : ncell) .');
IJK = [IJK{:}];

G.cells = struct('num',      numActiveCells,       ...
                 'numFaces', uint32(numCellFaces(active)), ...
                 'indexMap', uint32(find(active)),         ...
                 'ijkMap',   IJK(active,:));

I1 = find(cellFaces(:,3) > 0);
I2 = find(cellFaces(:,3) < 0);

neighbors                     = zeros([numFaces, 2]);
neighbors(cellFaces(I1,2), 1) = cellFaces(I1,1);
neighbors(cellFaces(I2,2), 2) = cellFaces(I2,1);

numNodes = uint8(accumarray(faceNodes(:,1), 1, [numFaces, 1]));
G.faces = struct('num',       numFaces, ...
                 'numNodes',  numNodes, ...
                 'neighbors', uint32(neighbors));

G.cellFaces =  int32(cellFaces(:,2:3));
G.faceNodes = uint32(faceNodes(:,2:3));
G.cartDims  = grdecl.cartDims(:) .';

%--------------------------------------------------------------------------
% Private helpers follow
%--------------------------------------------------------------------------


function [C, ncell, nnod] = cellCoordinates(grdecl)
nx = grdecl.cartDims(1);
ny = grdecl.cartDims(2);
nz = grdecl.cartDims(3);

nodeNo = reshape(1 : prod(2 .* [nx, ny, nz]), ...
                 2 .* [nx, ny, nz]);

n1 = nodeNo(1:2:2*nx, 1:2:2*ny, 1:2:2*nz);
n2 = nodeNo(2:2:2*nx, 1:2:2*ny, 1:2:2*nz);
n3 = nodeNo(1:2:2*nx, 2:2:2*ny, 1:2:2*nz);
n4 = nodeNo(2:2:2*nx, 2:2:2*ny, 1:2:2*nz);
n5 = nodeNo(1:2:2*nx, 1:2:2*ny, 2:2:2*nz);
n6 = nodeNo(2:2:2*nx, 1:2:2*ny, 2:2:2*nz);
n7 = nodeNo(1:2:2*nx, 2:2:2*ny, 2:2:2*nz);
n8 = nodeNo(2:2:2*nx, 2:2:2*ny, 2:2:2*nz);

nodeIndex = [n1(:), n2(:), n3(:), n4(:), ...
             n5(:), n6(:), n7(:), n8(:)];

pillarIx = reshape(1 : (nx+1)*(ny+1), nx + 1, ny + 1);
p1 = pillarIx(1:nx,   1:ny);
p2 = pillarIx(2:nx+1, 1:ny);
p3 = pillarIx(1:nx,   2:ny+1);
p4 = pillarIx(2:nx+1, 2:ny+1);

pCoord = repmat([p1(:), p2(:), p3(:), p4(:)], [nz, 2]);
pCoord = grdecl.COORD(pCoord(:), :);

zCoords = grdecl.ZCORN(nodeIndex(:));
propFac = (zCoords - pCoord(:,3)) ./ (pCoord(:,6) - pCoord(:,3));

xCoords = pCoord(:,1) + propFac.*(pCoord(:,4) - pCoord(:,1));
yCoords = pCoord(:,2) + propFac.*(pCoord(:,5) - pCoord(:,2));

C = [xCoords, yCoords, zCoords];
ncell = nx * ny * nz;
nnod  = 8 * ncell;


%--------------------------------------------------------------------------


function [F, activeF] = buildFaces(C, ncell, nnod)

N = reshape(1 : nnod, ncell, 8);
F = reshape([N(:,1), N(:,3), N(:,7), N(:,5), ...
             N(:,2), N(:,4), N(:,8), N(:,6), ...
             N(:,1), N(:,5), N(:,6), N(:,2), ...
             N(:,3), N(:,7), N(:,8), N(:,4), ...
             N(:,1), N(:,2), N(:,4), N(:,3), ...
             N(:,5), N(:,6), N(:,8), N(:,7)] .', 4, []) .';

T1 = F(:, [1, 2, 3]);
T2 = F(:, [3, 4, 1]);

normalT1 = cross(C(T1(:,1),:) - C(T1(:,2),:), ...
                 C(T1(:,3),:) - C(T1(:,2),:)) ./ 2;
normalT2 = cross(C(T2(:,1),:) - C(T2(:,2),:), ...
                 C(T2(:,3),:) - C(T2(:,2),:)) ./ 2;

areaT1 = sqrt(sum(normalT1 .^ 2, 2));
areaT2 = sqrt(sum(normalT2 .^ 2, 2));

activeF = areaT1 + areaT2 > sqrt(eps);

