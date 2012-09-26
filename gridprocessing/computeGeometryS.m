function G = computeGeometryS(G, varargin)
%Compute geometry of grid using single precision.
%
% SYNOPSIS:
%   G = computeGeometryS(G)
%   G = computeGeometryS(G, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - SAMSIM grid structure as described by grid_structure.
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display verbose output as
%                            the process progresses.  Possible values are
%                            TRUE and FALSE.  Default value is FALSE.
%
% RETURNS:
%   G - Grid structure with added fields:
%         - cells
%             - volumes
%             - centroids
%         - faces
%             - areas
%             - normals                !!! AREA WEIGHTED !!!
%             - centroids
%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas

%{
#COPYRIGHT#
%}

% $Id: computeGeometryS.m 1949 2009-03-31 10:16:08Z bska $

%% Setup
opt = struct('Verbose', false);
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;

numC = G.cells.num;
numF = G.faces.num;

faceNo  = rldecode(1:G.faces.num, double(G.faces.numNodes), 2) .';
p       = cumsum([0; double(G.faces.numNodes)]);
next    = (2:size(G.faceNodes, 1)+1) .';
next(p(2:end)) = p(1:end-1)+1; 
%%
% Divide each face into sub-triangles all having one node as
% pCenter = sum(nodes) /#nodes. Compute area-weighted normals,
% and add to obtain approx face-normals. Compute resulting areas
% and centroids.
if verbose,
   fprintf('Computing face normals/areas/centroids ... ');
   tic
end

llE = length(G.faceNodes);
%localEdge2Face = sparse(1:llE, double(G.faceNodes(:,1)), 1, llE, numF);
localEdge2Face=sparse(1:llE, faceNo, 1, llE, numF);

pCenters     = (localEdge2Face' * G.nodes.coords(G.faceNodes,:)) ...
                ./ double(G.faces.numNodes(:, [1, 1, 1]));
pCenters     = localEdge2Face*pCenters;
subNormals   = single(cross(G.nodes.coords(G.faceNodes(next),:) - ...
                            G.nodes.coords(G.faceNodes,:), ...
                            pCenters - G.nodes.coords(G.faceNodes,:)) ./ 2);
subAreas     = single(sqrt(sum(double(subNormals) .^ 2, 2)));
subCentroids = single((G.nodes.coords(G.faceNodes,:) + ...
                       G.nodes.coords(G.faceNodes(next),:) + pCenters) ./ 3);
clear pCenters llE;

faceNormals    = localEdge2Face' * double(subNormals);
subNormalSigns = int8(sign(sum(double(subNormals) .* (localEdge2Face * faceNormals), 2)));
faceAreas      = localEdge2Face.' * subAreas;
faceCentroids  = (localEdge2Face' * (repmat(double(subAreas) .* double(subNormalSigns), [1, 3]) .* ...
                     double(subCentroids))) ./ faceAreas(:, [1, 1, 1]);
clear subAreas;

tocif(verbose)

%%
% Divide each cell into sub-tetrahedra according to sub-triangles above,
% all having one node as cCenter = sum(faceCentroids) / #faceCentroids.

if verbose,
   fprintf('Computing cell volumes/centroids ... ');
   tic
end

cellVolumes   = zeros([numC, 1],'double');
cellCentroids = zeros([numC, 3],'double');

lastInx = uint32(0);
for c = 1 : numC,
   nF  = uint32(G.cells.numFaces(c));
   inx = (1:nF) + lastInx;

   faces        = G.cellFaces(inx,1);
   [triE, triF] = find(localEdge2Face(:,faces));
   nT           = length(triE);

   fCentroids = faceCentroids(faces,:);
   cCenter    = sum(fCentroids) ./ double(nF);

   relSubC     = double(subCentroids(triE,:)) - cCenter(ones([nT, 1]), :);
   %orientation = G.cellFaces(inx,2);
   orientation = 2*double(G.faces.neighbors(G.cellFaces(inx,1), 1)==c)-1;

   outNormals = double(subNormals    (triE, :)        ) .* ...
                double(subNormalSigns(triE, [1, 1, 1])) .* ...
                double(orientation   (triF, [1, 1, 1]));

   tVolumes   = abs((1/3) * sum(relSubC .* outNormals, 2));
   tCentroids = (3/4) * relSubC;

   volume      = sum(tVolumes);
   relCentroid = (tVolumes' * tCentroids) ./ volume;
   centroid    = relCentroid + cCenter;

   cellVolumes(c)     = double(volume);
   cellCentroids(c,:) = double(centroid);

   lastInx = lastInx + nF;
end
tocif(verbose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update grid
G.faces.areas     = faceAreas;     clear faceAreas;
G.faces.normals   = faceNormals;   clear faceNormals;
G.faces.centroids = faceCentroids; clear faceCentroids;

G.cells.volumes   = cellVolumes;   clear cellVolumes;
G.cells.centroids = cellCentroids; clear cellCentroids;
