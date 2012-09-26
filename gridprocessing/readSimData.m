function G = readSimData(caseName, varargin)
%Read SAM simulator grid input files.
%
% SYNOPSIS:
%   G = readSimData(caseName)
%   G = readSimData(caseName, maxIJK)
%
% PARAMETERS:
%   caseName - Simulation case base name.  String assumed to contain the
%              (relative) base name of a simulation case.
%
%              The grid of a simulation case is defined by three input
%              files:
%                 Geometry: [caseName, '-geom.dat']
%                 Topology: [caseName, '-topo.dat']
%                 IJK-map:  [caseName, '-map.dat']
%   maxIJK   - Vector, length 3, of maximal extents of cell logical,
%              cartesian indices.  OPTIONAL.
%              If not present (or empty), the cell IJK map is not read
%              and the G.cell.indexMap field will be empty.
%
% RETURNS:
%   G - Grid structure as detailed in grid_structure.
%
% EXAMPLE:
%   Assume that the input files for case 'example' are stored in
%   directory 'data' in the parent directory of the current working
%   directory (CD).  Assume furthermore that maximum extents of the
%   cell logical cartesian indices are 60, 220, and 85, respectively.
%   Then the case grid may be input as
%
%      caseName = ['..', filesep, 'data', filesep, 'example'];
%      maxIJK   = [60, 220, 85];
%      G = readSimData(caseName, maxIJK);
%
% SEE ALSO:
%   cartGrid, filesep, fopen, fscanf, textscan.

%{
#COPYRIGHT#
%}

% $Id: readSimData.m 1949 2009-03-31 10:16:08Z bska $

error(nargchk(1, 2, nargin));

%-------------------------------------------------------------
%% Read grid topology specification --------------------------
%
% File format is
%   [header -- inconsequential in MATLAB]
%   [blank]
%   [numCells numCellFaces numFaces numNodes]
%   [blank]
%   [numCells lines of
%       nLocFaces  locFace(1 : nLocFaces)], locFace indexed from zero.
%   [blank]
%   [numCellFaces lines of cellFace/orientation pairs]
%   [blank]
%   [numFaces lines of
%       nNodes globNode(1 : nNodes)], globNode indexed from zero.

loc = [caseName, '-topo.dat'];
[fid, msg] = fopen(loc, 'r');
if fid < 0, error(msg); end

header = fgetl(fid);

numCells     = fscanf(fid, '%d', 1);
numCellFaces = fscanf(fid, '%d', 1);
numFaces     = fscanf(fid, '%d', 1);
numNodes     = fscanf(fid, '%d', 1);

%% cellFaces

cellFaces       = zeros([numCellFaces, 3]);
numCellFacesVec = zeros([numCells, 1]);
for cellIx = 1 : numCells,
   numLocFaces = fscanf(fid, '%f', 1);
   numCellFacesVec(cellIx) = numLocFaces;

   cellFIx = fscanf(fid, '%f', numLocFaces) + 1;
   cellFaces(cellFIx, 1) = cellIx(ones([numLocFaces, 1]));
end

locF2F = fscanf(fid, '%f', [2, numCellFaces]) .';

orientations  = [-1, 1];
cellFaces(:,2) = locF2F(:,1) + 1;
cellFaces(:,3) = orientations(locF2F(:,2) + 1);

%% faceNodes
faceNodes       = zeros([6*numFaces, 3]);
numFaceEdgesVec = zeros([1*numFaces, 1]);
count = 0;
for faceIx = 1 : numFaces,
   numFNod = fscanf(fid, '%f', 1);
   nodeIx  = fscanf(fid, '%f', numFNod).' + 1;
   numFaceEdgesVec(faceIx) = numFNod;

   range = count + (1 : numFNod).';
   count = count + numFNod;

   faceNodes(range, 1) = faceIx(ones([numFNod, 1]));
   faceNodes(range, 2) = nodeIx;
   faceNodes(range, 3) = nodeIx([2:end, 1].');
end
faceNodes(count+1 : end, :) = [];

fclose(fid);


%-------------------------------------------------------------
%% Read grid geometry specification --------------------------
%
% File format is
%   [header -- inconsequential in MATLAB]
%   [blank]
%   [numGlobalNodes]
%   [blank]
%   [numGlobalNodes lines of X Y Z coordinates]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of unit face normals]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of face centroids]
%   [blank]
%   [numGlobalFaces]
%   [numGlobalFaces lines of face areas]
%   [blank]
%   [numGlobalCells]
%   [numGlobalCells lines of cell centroids]
%   [blank]
%   [numGlobalCells]
%   [numGlobalCells lines of cell volumes]

loc = [caseName, '-geom.dat'];
[fid, msg] = fopen(loc, 'r');
if fid < 0, error(msg); end

header = fgetl(fid);

numNodes      = fscanf(fid, '%f', 1);
coords        = fscanf(fid, '%f', [3, numNodes]) .';

numFaces      = fscanf(fid, '%f', 1);
normals       = fscanf(fid, '%f', [3, numFaces]) .';

numFaces      = fscanf(fid, '%f', 1);
faceCentroids = fscanf(fid, '%f', [3, numFaces]) .';

numFaces      = fscanf(fid, '%f', 1);
areas         = fscanf(fid, '%f', [1, numFaces]) .';

numCells      = fscanf(fid, '%f', 1);
cellCentroids = fscanf(fid, '%f', [3, numCells]) .';

numCells      = fscanf(fid, '%f', 1);
volumes       = fscanf(fid, '%f', [1, numCells]) .';

fclose(fid);

%-------------------------------------------------------------
%% Read grid cell IJK map specification ----------------------
%
% File format is
%   [numCell]
%   [blank]
%   [numGlobalNodes lines of cell logical I J K coordinates]

%if nargin > 1 && isnumeric(varargin{1}) && ~isempty(varargin{1}),
   loc = [caseName, '-map.dat'];
   [fid, msg] = fopen(loc, 'r');
   if fid < 0, error(msg); end

   NxNyNz   = fscanf(fid, '%d', [1 3]);
   numCells = fscanf(fid, '%f', 1);
   IJK      = fscanf(fid, '%f', [3, numCells]) .';
   IJK      = IJK + 1;     % CASENAME-map.dat indexes from zero

%   indexMap = IJK(:,1) + ...
%              varargin{1}(1)*((IJK(:,2) - 1)+ ...
%                              varargin{1}(2)*(IJK(:,3) - 1));
   indexMap = IJK(:,1) + ...
              NxNyNz(1)*((IJK(:,2) - 1)+ ...
                              NxNyNz(2)*(IJK(:,3) - 1));

   fclose(fid)
%else
%   IJK      = [];
%   indexMap = [];
%end


%% Find face-based neighbors (potentially SLOW)
neighbors = zeros(numFaces, 2);
for k = 1:numCellFaces,
   v = cellFaces(k,:);
   cellIx = v(1); faceIx = v(2); ori = v(3);
   if ori > 0
      neighbors(faceIx, 1) = cellIx;
   else
      neighbors(faceIx, 2) = cellIx;
   end
end


%-------------------------------------------------------------
%% Assemble structure ----------------------------------------
G.cells.num       = numCells;
G.cells.numFaces  = numCellFacesVec;
G.cells.volumes   = volumes;
G.cells.centroids = cellCentroids;
G.cells.indexMap  = indexMap;
G.cells.ijkMap    = IJK;

G.faces.neighbors = neighbors;
G.faces.num       = numFaces;
G.faces.numNodes  = numFaceEdgesVec;
G.faces.areas     = areas;
G.faces.normals   = normals .* areas(:, [1, 1, 1]);
G.faces.centroids = faceCentroids;

G.nodes.num    = numNodes;
G.nodes.coords = coords;

%G.cellFaces = cellFaces(:,2:3);
G.cellFaces = cellFaces;
G.faceNodes = faceNodes(:,2:3);
G.cartDims = NxNyNz;
