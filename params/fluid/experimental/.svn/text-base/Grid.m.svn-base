classdef Grid < handle
    
    properties(Access=public, Abstract=true)
        type
        cells
        faces
        nodes
        cellFaces
        faceNodes
        cartDims
    end
    
    methods(Access=public)
        function self = computeGeometry(self, varargin)
            % COMPUTEGEOMETRY -- Compute geometry of grid for mimetic system
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

            % $Id: computeGeometry.m 464 2008-07-02 10:22:50Z jrn $

            %% Setup
            opt = struct('Verbose', false);
            opt = merge_options(opt, varargin{:});
            verbose = opt.Verbose;

            numC = self.cells.num;
            numF = self.faces.num;

            %%
            % Divide each face into sub-triangles all having one node as
            % pCenter = sum(nodes) /#nodes. Compute area-weighted normals,
            % and add to obtain approx face-normals. Compute resulting areas
            % and centroids.
            if verbose,
                fprintf('Computing face normals/areas/centroids ... ');
                tic
            end

            llE = length(self.faceNodes);
            %localEdge2Face = sparse(1:llE, double(self.faceNodes(:,1)), 1, llE, numF);
            localEdge2Face=sparse(1:llE, rldecode((1:self.faces.num)',double(self.faces.numNodes)), 1, llE, numF);

            pCenters     = (localEdge2Face' * self.nodes.coords(self.faceNodes(:,1),:)) ...
                ./ double(self.faces.numNodes(:, [1, 1, 1]));
            pCenters     = localEdge2Face*pCenters;
            subNormals   = single(cross(self.nodes.coords(self.faceNodes(:,2),:) - ...
                self.nodes.coords(self.faceNodes(:,1),:), ...
                pCenters - self.nodes.coords(self.faceNodes(:,1),:)) ./ 2);
            subAreas     = single(sqrt(sum(double(subNormals) .^ 2, 2)));
            subCentroids = single((self.nodes.coords(self.faceNodes(:,1),:) + ...
                self.nodes.coords(self.faceNodes(:,2),:) + pCenters) ./ 3);
            clear pCenters llE;

            faceNormals    = localEdge2Face' * double(subNormals);
            subNormalSigns = int8(sign(sum(double(subNormals) .* (localEdge2Face * faceNormals), 2)));
            faceAreas      = sqrt(sum(faceNormals .^ 2, 2));
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
                nF  = uint32(self.cells.numFaces(c));
                inx = (1:nF) + lastInx;

                faces        = self.cellFaces(inx,1);
                [triE, triF] = find(localEdge2Face(:,faces));
                nT           = length(triE);

                fCentroids = faceCentroids(faces,:);
                cCenter    = sum(fCentroids) ./ double(nF);

                relSubC     = double(subCentroids(triE,:)) - cCenter(ones([nT, 1]), :);
                orientation = self.cellFaces(inx,2);

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
            self.faces.areas     = faceAreas;     clear faceAreas;
            self.faces.normals   = faceNormals;   clear faceNormals;
            self.faces.centroids = faceCentroids; clear faceCentroids;

            self.cells.volumes   = cellVolumes;   clear cellVolumes;
            self.cells.centroids = cellCentroids; clear cellCentroids;

        end

    end
end
