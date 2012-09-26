function G = computeGeometry(G, varargin)
%Compute geometry of grid.
%
% SYNOPSIS:
%   G = computeGeometry(G)
%   G = computeGeometry(G, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
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

% $Id: computeGeometry.m 2982 2009-10-12 13:35:19Z jrn $

%% Setup
assert(size(G.faceNodes, 2)==1);
opt     = struct('verbose', false);
opt     = merge_options(opt, varargin{:});

numC    = G.cells.num;
numF    = G.faces.num;

if size(G.nodes.coords,2) == 3,
   faceNo  = rldecode(1:G.faces.num, double(G.faces.numNodes), 2) .';
   p       = cumsum([0; double(G.faces.numNodes)]);
   next    = (2:size(G.faceNodes, 1)+1) .';
   next(p(2:end)) = p(1:end-1)+1;
   %%
   % Divide each face into sub-triangles all having one node as
   % pCenter = sum(nodes) /#nodes. Compute area-weighted normals,
   % and add to obtain approx face-normals. Compute resulting areas
   % and centroids.
   dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
   ticif (opt.verbose)

   llE = length(G.faceNodes);
   localEdge2Face = sparse(1 : llE, faceNo, 1, llE, numF);

   pCenters     = bsxfun(@rdivide, ...
                         localEdge2Face.' * G.nodes.coords(G.faceNodes,:), ...
                         double(G.faces.numNodes));
   pCenters     = localEdge2Face * pCenters;
   subNormals   = cross(G.nodes.coords(G.faceNodes(next),:) - ...
                        G.nodes.coords(G.faceNodes,:), ...
                        pCenters - G.nodes.coords(G.faceNodes,:)) ./ 2;
   subAreas     = sqrt(sum(subNormals .^ 2, 2));
   subCentroids = (G.nodes.coords(G.faceNodes,:) + ...
                   G.nodes.coords(G.faceNodes(next),:) + pCenters) ./ 3;
   clear pCenters llE faceNo

   faceNormals    = localEdge2Face.' * subNormals;
   faceAreas      = localEdge2Face.' * subAreas;
   subNormalSigns = sign(sum(subNormals .* (localEdge2Face * faceNormals), 2));
   faceCentroids  = bsxfun(@rdivide,                                 ...
                           localEdge2Face.' * ...
                              bsxfun(@times, subAreas, subCentroids), ...
                           faceAreas);
   clear subAreas
   
   tocif(opt.verbose)

   %%
   % Divide each cell into sub-tetrahedra according to sub-triangles above,
   % all having one node as cCenter = sum(faceCentroids) / #faceCentroids.

   dispif(opt.verbose, 'Computing cell volumes and centroids...\t\t');
   ticif (opt.verbose);

   cellVolumes   = zeros([numC, 1]);
   cellCentroids = zeros([numC, 3]);

   lastInx = 0;
   for c = 1 : numC,
      nF  = double(G.cells.numFaces(c));
      inx = (1 : nF) + lastInx;

      faces        = G.cellFaces(inx,1);
      [triE, triF] = find(localEdge2Face(:,faces));

      fCentroids = faceCentroids(faces,:);
      cCenter    = sum(fCentroids) ./ double(nF);

      relSubC    = bsxfun(@minus, subCentroids(triE,:), cCenter);

      % The normal of a face f is directed from cell G.faces.neighbors(f,1) to
      % cell G.faces.neighbors(f,2).   If cell c is in the second column for
      % face f, then the nomal must be multiplied by -1 to be an outer normal.
      orientation = 2*double(G.faces.neighbors(G.cellFaces(inx,1), 1) == c)-1;

      outNormals = bsxfun(@times,             ...
                          subNormals(triE,:), ...
                          subNormalSigns(triE) .* orientation(triF));

      tVolumes   = abs((1/3) * sum(relSubC .* outNormals, 2));
      tCentroids = (3/4) * relSubC;

      volume      = sum(tVolumes);
      relCentroid = (tVolumes' * tCentroids) ./ volume;
      centroid    = relCentroid + cCenter;

      cellVolumes(c)     = volume;
      cellCentroids(c,:) = centroid;

      lastInx = lastInx + nF;
   end
   tocif(opt.verbose)

elseif size(G.nodes.coords,2)==2

  dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
  ticif (opt.verbose)
 
  faceEdges     = reshape(G.faceNodes,2,[])';
  edgeLength    = G.nodes.coords(faceEdges(:,2),:)-G.nodes.coords(faceEdges(:,1),:);
  faceAreas     = sqrt(sum(edgeLength.*edgeLength,2));
  faceCentroids = (G.nodes.coords(faceEdges(:,2),:)+ ...
                   G.nodes.coords(faceEdges(:,1),:))/2;
  faceNormals   = [edgeLength(:,2),-edgeLength(:,1)];

  tocif(opt.verbose)
  dispif(opt.verbose, 'Computing cell volumes and centroids ...\t\t');
  ticif (opt.verbose)
  
  cellVolumes   = zeros([numC, 1]);
  cellCentroids = zeros([numC, 2]);

  %cellno = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';

  lastInx = 0;
  for i=1:numC,
    nF  = double(G.cells.numFaces(i));
    inx = (1 : nF) + lastInx;
     
    faces        = G.cellFaces(inx,1);
    cCenter     = sum(faceCentroids(faces, :))/numel(faces);
    a           = bsxfun(@minus, G.nodes.coords(faceEdges(faces,2),:), cCenter);
    b           = bsxfun(@minus, G.nodes.coords(faceEdges(faces,1),:), cCenter);
    subArea     = a(:,1).*b(:,2)-a(:,2).*b(:,1);  
    subArea     = 0.5*abs(subArea);
    subCentroid = bsxfun(@plus, cCenter,  2*faceCentroids(faces,:))/3;
    cellVolumes(i) = sum(subArea); 
    cellCentroids(i, :) = sum(bsxfun(@times, subArea, subCentroid))/cellVolumes(i);
    
    lastInx = lastInx + nF;
  end

  tocif(opt.verbose)

else 
  assert(false);  
end

%% Update grid
G.faces.areas     = faceAreas;
G.faces.normals   = faceNormals;
G.faces.centroids = faceCentroids;

G.cells.volumes   = cellVolumes;
G.cells.centroids = cellCentroids;
