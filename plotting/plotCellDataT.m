function [ht,hp] = plotCellDataT(G, data, subDomain, varargin)
%Plot triangulated exterior grid faces, variable colouring, to current axes.
%
% SYNOPSIS:
%   plotCellDataT(G, data);
%   plotCellDataT(G, data, subDomain);
%   plotCellDataT(G, data, [], 'PropertyName', PropertyValue, ...);
%   plotCellDataT(G, data, subDomain, 'PropertyName', PropertyValue, ...);
%   [ht hp] = plotCellDataT(...)
%
% PARAMETERS:
%   G         - Grid-structure
%   data      - Vector of length G.cells.num
%   subDomain        - Vector of cell-indices defining sub-grid.
%                      OPTIONAL. Default value = [], which is equivalent to (1:G.cells.num)'.
%   patch properties - List of property names/property values. OPTIONAL.
%                      Default value = {'FaceColor', 'flat', 'EdgeColor', 'none'}
%
% RETURNS:
%   ht - handle to triangular patch-object of triangulated grid surface
%   hp - Structure of handles to topological unique patch-objects, i.e. hp(1) refers
%        to all triangles, hp(2) to all quads etc. It has the following fields:
%         - patch   :  handle to current patch object
%         - polyN   :  current polygon-type (integer)
%
% COMMENTS:
%   This plotting routine is better than plotCellData for strongly curved
%   faces.
%
% REFERENCE:
%   See also plotCellData, plotGrid, plotGridT

%{
#COPYRIGHT#
%}

% $Id: plotCellDataT.m 1952 2009-03-31 10:30:47Z bska $

if nargin < 3,
   subDomain=[];                           % actually subDomain = (1:G.cells.num)';
end

if ~isnumeric(subDomain),
   error('Given subdomain not numeric' )
end

%% Put properties/values in structure patchPropEdge and patchPropFace
patchPropFace = struct('FaceColor', 'flat', 'EdgeColor', 'none');
patchPropEdge = struct('FaceColor', 'none', 'EdgeColor', 'none');

if ~isempty(varargin),
   addProps = struct(varargin{:});
   fNames = fieldnames(addProps);
   for k = 1:length(fNames)
      fName = fNames{k};
      if ~strcmp(fName,'EdgeColor')
         patchPropFace.(fName) = addProps.(fName);
      end
      if ~strcmp(fName,'FaceColor')
         patchPropEdge.(fName) = addProps.(fName);
      end
   end
end

%% Find grid/subgrid boundary faces
V = G.faces.neighbors;
if ~isempty(subDomain),
   ind = zeros(G.cells.num,1);
   ind(subDomain) = subDomain;
   V = G.faces.neighbors;
   pos = find(V);
   V(pos) = ind(V(pos));
   bndryIndex = and( V(:,1).*V(:,2) == 0, V(:,1)+V(:,2) ~= 0) ;
else
   bndryIndex = ( (V(:,1).*V(:,2) ) == 0 );
end

bndryFaces = find(bndryIndex);
data = data(:);
bndryData = data( max(V(bndryFaces,:), [], 2) );
clear V;

faceNo = rldecode((1:G.faces.num)',double(G.faces.numNodes));
localBndryIndex = find(bndryIndex(faceNo));
localBndryEdges = [faceNo(localBndryIndex), G.faceNodes(localBndryIndex,:)];
clear faceNo localBndryIndex;

%% Re-index
ind = zeros(size(bndryIndex));
ind(bndryFaces)=(1:length(bndryFaces))';
faceNodes = [ind(localBndryEdges(:,1)) localBndryEdges(:,2:3)];
clear localBndryEdges bndryIndex ind;

%% Sparse mapping from exterior edge to face
C = sparse((1 : size(faceNodes,1)).', double(faceNodes(:,1)), 1, ...
           size(faceNodes,1), length(bndryFaces));

%% Compute face "centers"
centers = (C'*G.nodes.coords(faceNodes(:,2),:))./(sum(C)'*[1 1 1]);

%% plot Triangel faces
numN = G.nodes.num;
ht = patch('Faces', [faceNodes(:,2:3) faceNodes(:,1)+numN], ...
           'Vertices',[G.nodes.coords; centers], ...
           'FaceVertexCData', bndryData(faceNodes(:,1)), ...
           patchPropFace);

%% plot Ploygon edges
numNodes = double(G.faces.numNodes(bndryFaces));
instances = unique(numNodes);                % occuring polygons

hp = struct;
for k = 1 : length(instances),
   n = instances(k);
   hp(k).polyN = n;
   ind = find(numNodes == n);
   nodes = faceNodes(sum(C(:,ind), 2)~=0, 2); %#ok
   hp(k).patch = patch('Faces', reshape(nodes,n,length(nodes)/n)', ...
                       'Vertices', G.nodes.coords, ...
                       patchPropEdge);
end

set(gca,'ZDir','reverse')
