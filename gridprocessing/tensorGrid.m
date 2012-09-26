function G = tensorGrid(dx, dy, varargin)
%Construct Cartesian grid with variable physical cell sizes.
%
% SYNOPSIS:
%   G = tensorGrid(dx, dy)
%   G = tensorGrid(dx, dy, 'depthz, z)
%   G = tensorGrid(dx, dy, dz)
%   G = tensorGrid(dx, dy, dz, 'depthz', z)
%
% PARAMETERS:
%   dx,dy,dz - Vectors giving cell sizes, in units of meters, of individual
%              coordinate directions.  Specifically, the grid cell at
%              logical location (I,J,K) will have a physical dimension of
%              [dx(I), dy(J), dz(K)] (meters).
%
%   depthz   - Depth, in units of meters, at which upper reservoir nodes
%              are encountered.  Assumed to be a
%              (NUMEL(dx)+1)-by-(NUMEL(dy)+1) array of nodal depths.
%
%              OPTIONAL.
%              Default value: depthz = ZEROS([numel(dx)+1, numel(dy)+1])
%                             (i.e., top of reservoir at zero depth).
%
% RETURNS:
%   G - Grid structure mostly as detailed in grid_structure, though lacking the
%       fields
%         - G.cells.volumes
%         - G.cells.centroids
%
%         - G.faces.areas
%         - G.faces.normals
%         - G.faces.centroids
%
%       These fields may be computed using the function computeGeometry.
%
%       There is, however, an additional field not described in grid_structure:
%
%         - cartDims -- A length 2 or 3 vector giving number of cells in each
%                       coordinate direction.  In other words
%
%                         cartDims == [NUMEL(dx), NUMEL(dy), NUMEL(dz)].
%
%  New:  G.cellFaces(:,2) contains integers 1-6 corresponding to directions
%        W, E, S, N, T, B.
%
% SEE ALSO:
%   grid_structure, computeGeometry

%{
#COPYRIGHT#
%}

% $Date: 2009-10-14 08:03:36 +0200 (on, 14 okt 2009) $
% $Revision: 2993 $

   
   
   
   if mod(nargin, 2), %3D
      G = tensorGrid3D(dx, dy, varargin{:});
   else
      G = tensorGrid2D(dx, dy, varargin{:});
   end
   
end

function G = tensorGrid3D(dx, dy, dz, varargin)
   error(nargchk(3, 5, nargin, 'struct'));
   
   sx = numel(dx); sy = numel(dy); sz = numel(dz);
   celldim = [sx, sy, sz];
   
   opt = struct('depthz', zeros(sx+1, sy+1));
   opt = merge_options(opt, varargin{:});
   if numel(opt.depthz) ~= (sx+1) * (sy+1),
      error(msgid('DepthZ:WrongSize'), ...
         'Input argument ''depthz'' is wrongly sized.')
   end
   
   
   
   
   numC = sx * sy * sz;             % Number of cells.
   numN = (sx+1) * (sy+1) * (sz+1); % Number of nodes.
   
   numFX = (sx+1) * sy * sz;        % Number of faces parallel to yz-plane.
   numFY = sx * (sy+1) * sz;        % Number of faces parallel to xz-plane.
   numFZ = sx * sy * (sz+1);        % Number of faces parallel to xy-plane.
   numF  = numFX + numFY + numFZ;
   
   %--------------------------------------------------------------------------
   % Nodes/Coordinates -------------------------------------------------------
   
   x = [0; cumsum(dx(:))];
   y = [0; cumsum(dy(:))];
   z = [0; cumsum(dz(:))];
   
   [xCoord, yCoord, zCoord] = ndgrid(x, y, z);
   zCoord = zCoord + repmat(reshape(opt.depthz, sx+1, sy+1), [1, 1, sz+1]);
   
   coords = [xCoord(:), yCoord(:), zCoord(:)];
   
   %--------------------------------------------------------------------------
   %% Generate face-edges ----------------------------------------------------
   
   % Node index matrix
   N = reshape(1 : numN, [sx+1, sy+1, sz+1]);
   
   %--------------------------------------------------------------------------
   % x-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx+1, 1:sy  , 1:sz  ), 1, []);
   NF2 = reshape(N(1:sx+1, 2:sy+1, 1:sz  ), 1, []);
   NF3 = reshape(N(1:sx+1, 2:sy+1, 2:sz+1), 1, []);
   NF4 = reshape(N(1:sx+1, 1:sy  , 2:sz+1), 1, []);
   
   faceNodesX = reshape([NF1; NF2; NF3; NF4], [], 1);
   
   %--------------------------------------------------------------------------
   % y-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy+1, 1:sz  ), 1, []);
   NF2 = reshape(N(1:sx  , 1:sy+1, 2:sz+1), 1, []);
   NF3 = reshape(N(2:sx+1, 1:sy+1, 2:sz+1), 1, []);
   NF4 = reshape(N(2:sx+1, 1:sy+1, 1:sz  ), 1, []);
   
   faceNodesY = reshape([NF1; NF2; NF3; NF4], [], 1);
   
   %--------------------------------------------------------------------------
   % z-faces -----------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy,   1:sz+1), 1, []);
   NF2 = reshape(N(2:sx+1, 1:sy,   1:sz+1), 1, []);
   NF3 = reshape(N(2:sx+1, 2:sy+1, 1:sz+1), 1, []);
   NF4 = reshape(N(1:sx  , 2:sy+1, 1:sz+1), 1, []);
   
   faceNodesZ = reshape([NF1; NF2; NF3; NF4], [], 1);
   
   %--------------------------------------------------------------------------
   % Assemble grid_structure faceNodes structure -----------------------------
   %
   faceNodes = [faceNodesX; ...
      faceNodesY; ...
      faceNodesZ];
   
   clear -regexp ^faceNodes. ^N ^NF.
   %--------------------------------------------------------------------------
   %% Generate cell-faces ----------------------------------------------------
   
   foffset = 0;
   % Face index matrices
   FX = reshape(foffset + (1:numFX), sx+1, sy  , sz  ); foffset = foffset + numFX;
   FY = reshape(foffset + (1:numFY), sx  , sy+1, sz  ); foffset = foffset + numFY;
   FZ = reshape(foffset + (1:numFZ), sx  , sy  , sz+1);
   
   F1 = reshape(FX(1:sx  , :, :), 1, []); %W == 1
   F2 = reshape(FX(2:sx+1, :, :), 1, []); %E == 2
   
   F3 = reshape(FY(:, 1:sy  , :), 1, []); %S == 3
   F4 = reshape(FY(:, 2:sy+1, :), 1, []); %N == 4
   
   F5 = reshape(FZ(:, :, 1:sz  ), 1, []); %T == 5
   F6 = reshape(FZ(:, :, 2:sz+1), 1, []); %B == 6
   
   cellFaces = [reshape([F1; F2; F3; F4; F5; F6], [], 1), ...
      kron(ones([numC, 1]), [ 1, 2, 3, 4, 5, 6]')];
   
   clear -regexp ^F.
   
   %--------------------------------------------------------------------------
   %% Generate neighbors -----------------------------------------------------
   
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
   
   clear -regexp ^N..
   
   %--------------------------------------------------------------------------
   %% Assemble structure -----------------------------------------------------
   
   G.type  = 'tensorgrid';
   G.cells = struct('num',      numC,                         ...
      'numFaces', repmat(uint8(6), [numC, 1]),  ...
      'facePos',  uint32(1:6:(numC+1)*6)',      ...
      'indexMap', uint32(1 : numC)');
   
   G.faces = struct('num',       numF,                        ...
      'numNodes',  repmat(uint8(4), [numF, 1]), ...
      'nodePos',   uint32(1:4:(numF+1)*4)',     ...
      'neighbors', uint32(neighbors),           ...
      'tag',       zeros(numF, 1, 'uint8'));
   
   G.nodes = struct('num', uint32(numN), 'coords', coords);
   
   G.cellFaces = uint32(cellFaces);
   G.faceNodes = uint32(faceNodes);
   G.cartDims  = celldim;
end


function G =tensorGrid2D(dx, dy, varargin)
   error(nargchk(2, 4, nargin, 'struct'));
   
   sx = numel(dx); sy = numel(dy);
   celldim = [sx, sy];
   
   numC = sx * sy;             % Number of cells.
   numN = (sx+1) * (sy+1)   ; % Number of nodes.

   numFX = (sx+1) * sy;        % Number of faces parallel to y-axis.
   numFY = sx * (sy+1);        % Number of faces parallel to x-axis.
   numF  = numFX + numFY;

   %-----------------------------------------------------------------------
   % Nodes/Coordinates ----------------------------------------------------
   x = [0; cumsum(dx(:))];
   y = [0; cumsum(dy(:))];

   [xCoord, yCoord] = ndgrid(x, y);
   coords = [xCoord(:), yCoord(:)];

   %-----------------------------------------------------------------------
   %% Generate face-edges -------------------------------------------------

   % Node index matrix
   N = reshape(1 : numN, [sx+1, sy+1]);

   %-----------------------------------------------------------------------
   % x-faces --------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx+1, 1:sy  ), 1, []);
   NF2 = reshape(N(1:sx+1, 2:sy+1), 1, []);

   faceNodesX = reshape([NF1; NF2], [], 1);

   %-----------------------------------------------------------------------
   % y-faces --------------------------------------------------------------
   %
   NF1 = reshape(N(1:sx  , 1:sy+1), 1, []);
   NF2 = reshape(N(2:sx+1, 1:sy+1), 1, []);

   faceNodesY = reshape([NF2; NF1], [], 1);
   % Note:
   % Nodes need to be reversed to obtain normals
   % pointing in positive i-direction in computeGeometry.

   %-----------------------------------------------------------------------
   % Assemble grid_structure faceNodes structure --------------------------
   %
   faceNodes = [faceNodesX; ...
                faceNodesY];

   clear -regexp ^faceNodes. ^N ^NF.
   %-----------------------------------------------------------------------
   %% Generate cell-faces -------------------------------------------------

   foffset = 0;
   % Face index matrices
   FX = reshape(foffset + (1:numFX), sx+1, sy  ); foffset = foffset + numFX;
   FY = reshape(foffset + (1:numFY), sx  , sy+1);

   F1 = reshape(FX(1:sx  , :), 1, []); %W == 1
   F2 = reshape(FX(2:sx+1, :), 1, []); %E == 2

   F3 = reshape(FY(:, 1:sy  ), 1, []); %S == 3
   F4 = reshape(FY(:, 2:sy+1), 1, []); %N == 4

   cellFaces = [reshape([F1; F2; F3; F4], [], 1), ...
                kron(ones([numC, 1]), [ 1, 2, 3, 4]')];

   clear -regexp ^F.

   %-----------------------------------------------------------------------
   %% Generate neighbors --------------------------------------------------

   % Cell index matrix
   C = zeros([sx+2, sy+2]);
   C(2:sx+1, 2:sy+1) = reshape(1:numC, [sx, sy]);

   NX1 = reshape(C( 1:sx+1, 2:sy+1), [], 1);
   NX2 = reshape(C( 2:sx+2, 2:sy+1), [], 1);

   NY1 = reshape(C( 2:sx+1, 1:sy+1), [], 1);
   NY2 = reshape(C( 2:sx+1, 2:sy+2), [], 1);

   neighbors = [[NX1, NX2]; [NY1, NY2]];

   clear -regexp ^N..

   %-----------------------------------------------------------------------
   %% Assemble structure --------------------------------------------------

   G.type  = 'tensorgrid';
   G.cells = struct('num',      numC,                         ...
                    'numFaces', repmat(uint8(4), [numC, 1]),  ...
                    'facePos',  uint32(1:4:(numC+1)*4)',      ...
                    'indexMap', uint32(1 : numC)');

   G.faces = struct('num',       numF,                        ...
                    'numNodes',  repmat(uint8(2), [numF, 1]), ...
                    'nodePos',   uint32(1:2:(numF+1)*2)',     ...
                    'neighbors', uint32(neighbors),           ...
                    'tag',       zeros(numF, 1, 'uint8'));

   G.nodes = struct('num', uint32(numN), 'coords', coords);

   G.cellFaces = uint32(cellFaces);
   G.faceNodes = uint32(faceNodes);
   G.cartDims  = celldim;
end

