
classdef TensorGrid < Grid
    properties(Access=public)
        type
        cells
        faces
        nodes
        cellFaces
        faceNodes
        cartDims


    end
    methods(Access=public)
        function self = TensorGrid(dx, dy, dz, varargin)
            % TENSORGRID -- Construct cartesian grid with variable
            %               dx, dy, dz in physical space
            %
            % SYNOPSIS:
            %   G = tensorGrid(dx, dy, dz, depthz);
            %
            % PARAMETERS:
            %   dx,dy,dz  -
            %
            % RETURNS:
            %   G - Grid structure mostly as detailed in grid_structure, though lacking the
            %       fields
            %         - G.cells.volumes
            %         - G.cells.centroids
            %         - G.cells.ijkMap
            %
            %         - G.faces.areas
            %         - G.faces.normals
            %         - G.faces.centroids
            %
            %       These fields may be computed using the function computeGeometry.
            %
            %       There is, however, an additional field not described in grid_structure:
            %
            %         - cartDims -- A length 3 vector giving number of cells in each
            %                       coordinate direction.  In other words
            %
            %                               cartDims == celldim .
            %
            % SEE ALSO:
            %   grid_structure, computeGeometry

            % $Id: tensorGrid.m 213 2008-05-13 17:48:19Z knl $

            %%
            error(nargchk(3, 4, nargin, 'struct'));

            sx = numel(dx); sy = numel(dy); sz = numel(dz);
            celldim = [sx, sy, sz];

            if nargin == 3,
                depthz = zeros(sx+1, sy+1);
            else
                depthz = varargin{1};
            end

            numC = sx*sy*sz;                            % number of cells
            numN = (sx+1)*(sy+1)*(sz+1);                % number of nodes

            numFX = (sx+1)*sy*sz;                       % number of faces parallel to yz-plane
            numFY = sx*(sy+1)*sz;                       % number of faces parallel to xz-plane
            numFZ = sx*sy*(sz+1);                       % number of faces parallel to xy-plane
            numF = numFX + numFY + numFZ;

            %------------------------------------------------------------
            % Nodes/Coordinates -----------------------------------------

            x = [0; cumsum(dx(:))];
            y = [0; cumsum(dy(:))];
            z = [0; cumsum(dz(:))];

            [xCoord, yCoord, zCoord] = ndgrid(x, y, z);
            zCoord = zCoord + repmat(reshape(depthz, sx+1, sy+1), [1,1,sz+1]);

            coords = [xCoord(:), yCoord(:), zCoord(:)];

            %-------------------------------------------------------------
            %% Generate face-edges ---------------------------------------

            % Node index matrix
            N = reshape(1 : numN, [sx+1, sy+1, sz+1]);

            %-------------------------------------------------------------
            % x-faces ----------------------------------------------------
            %
            NF1 = reshape(N(1:sx+1, 1:sy  , 1:sz  ), 1, []);
            NF2 = reshape(N(1:sx+1, 2:sy+1, 1:sz  ), 1, []);
            NF3 = reshape(N(1:sx+1, 2:sy+1, 2:sz+1), 1, []);
            NF4 = reshape(N(1:sx+1, 1:sy  , 2:sz+1), 1, []);

            foffset = 0;
            faceNodesX = [kron(foffset + (1 : numFX).', ones([4, 1])), ...
                reshape([NF1; NF2; NF3; NF4], [], 1),        ...
                reshape([NF2; NF3; NF4; NF1], [], 1)];

            %-------------------------------------------------------------
            % y-faces ----------------------------------------------------
            %
            NF1 = reshape(N(1:sx  , 1:sy+1, 1:sz  ), 1, []);
            NF2 = reshape(N(1:sx  , 1:sy+1, 2:sz+1), 1, []);
            NF3 = reshape(N(2:sx+1, 1:sy+1, 2:sz+1), 1, []);
            NF4 = reshape(N(2:sx+1, 1:sy+1, 1:sz  ), 1, []);

            foffset = foffset + numFX;
            faceNodesY = [kron(foffset + (1 : numFY).', ones([4, 1])), ...
                reshape([NF1; NF2; NF3; NF4], [], 1),        ...
                reshape([NF2; NF3; NF4; NF1], [], 1)];

            %-------------------------------------------------------------
            % z-faces ----------------------------------------------------
            %
            NF1 = reshape(N(1:sx  , 1:sy,   1:sz+1), 1, []);
            NF2 = reshape(N(2:sx+1, 1:sy,   1:sz+1), 1, []);
            NF3 = reshape(N(2:sx+1, 2:sy+1, 1:sz+1), 1, []);
            NF4 = reshape(N(1:sx  , 2:sy+1, 1:sz+1), 1, []);

            foffset = foffset + numFY;
            faceNodesZ = [kron(foffset + (1 : numFZ).', ones([4, 1])), ...
                reshape([NF1; NF2; NF3; NF4], [], 1),        ...
                reshape([NF2; NF3; NF4; NF1], [], 1)];

            %-------------------------------------------------------------
            % Assemble grid_structure faceNodes structure --------------------
            %
            faceNodes = [faceNodesX; ...
                faceNodesY; ...
                faceNodesZ];

            clear -regexp ^faceNodes. ^N ^NF.
            %-------------------------------------------------------------
            %% Generate cell-faces --------------------------------------

            foffset = 0;
            % Face index matrices
            FX = reshape(foffset + (1:numFX), sx+1, sy  , sz  ); foffset = foffset + numFX;
            FY = reshape(foffset + (1:numFY), sx  , sy+1, sz  ); foffset = foffset + numFY;
            FZ = reshape(foffset + (1:numFZ), sx  , sy  , sz+1);

            F1 = reshape(FX(1:sx  , :, :), 1, []); %E
            F2 = reshape(FX(2:sx+1, :, :), 1, []); %W

            F3 = reshape(FY(:, 1:sy  , :), 1, []); %S
            F4 = reshape(FY(:, 2:sy+1, :), 1, []); %N

            F5 = reshape(FZ(:, :, 1:sz  ), 1, []); %T
            F6 = reshape(FZ(:, :, 2:sz+1), 1, []); %B

            cellFaces = [kron((1 : numC).', ones([6, 1])),         ...
                reshape([F1; F2; F3; F4; F5; F6], [], 1), ...
                kron(ones([numC, 1]), [-1, 1, -1, 1, -1, 1]')];

            clear -regexp ^F.

            %-------------------------------------------------------------
            %% Generate neighbors ----------------------------------------

            % Cell index matrix
            C = zeros([sx+2, sy+2, sz+2]);
            C(2:sx+1, 2:sy+1, 2:sz+1) = reshape(1:numC, [sx, sy, sz]);

            NX1 = reshape(C( 1:sx+1, 2:sy+1, 2:sz+1), [], 1);
            NX2 = reshape(C( 2:sx+2, 2:sy+1, 2:sz+1), [], 1);

            NY1 = reshape(C( 2:sx+1, 1:sy+1, 2:sz+1), [], 1);
            NY2 = reshape(C( 2:sx+1, 2:sy+2, 2:sz+1), [], 1);

            NZ1 = reshape(C( 2:sx+1, 2:sy+1, 1:sz+1), [], 1);
            NZ2 = reshape(C( 2:sx+1, 2:sy+1, 2:sz+2), [], 1);

            neighbors = [ [NX1, NX2]; ...
                [NY1, NY2]; ...
                [NZ1, NZ2] ];

            nE= (1   :sx+1:numFX)';
            nW= (sx+1:sx+1:numFX)';

            n = reshape(1:numFY, sx, sy+1, sz);
            nS= n(:, 1 ,:); nS=nS(:)+numFX;
            nN= n(:,end,:); nN=nN(:)+numFX;

            n = reshape(1:numFZ, sx, sy, sz+1);
            nT= n(:,:, 1 ); nT=nT(:)+numFX+numFY;
            nB= n(:,:,end); nB=nB(:)+numFX+numFY;

            tag = zeros(numF, 1);
            tag(nE) = 1;
            tag(nW) = 2;
            tag(nS) = 3;
            tag(nN) = 4;
            tag(nT) = 5;
            tag(nB) = 6;

            clear -regexp ^N..
            clear -regexp ^n.$

            %-------------------------------------------------------------
            %% Assemble structure ----------------------------------------

            self.type  = 'tensorgrid';
            self.cells = struct('num',      numC,                         ...
                'numFaces', repmat(uint8(6), [numC, 1]),  ...
                'indexMap', uint32(1 : numC)');

            self.faces = struct('num',       numF,                        ...
                'numNodes',  repmat(uint8(4), [numF, 1]), ...
                'neighbors', uint32(neighbors),           ...
                'tag',       tag);

            self.nodes = struct('num', uint32(numN), 'coords', coords);

            self.cellFaces =  int32(cellFaces(:,2:3));
            self.faceNodes = uint32(faceNodes(:,2:3));
            self.cartDims  = celldim;
        end
    end
end