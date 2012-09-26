function W = addWell000(G, rock, W, varargin)
error('Function ''putWell000'' is obsolete')
%{
%function W = addWell(G, rock, W, cellInx, type, val, radius, direction)
% 
%
% PUTWELL -- PUT WELL, SIMPLE SOURCE OR BOUNDARY CONDITIONS INTO
% WELL-STRUCTURE
%
% SYNOPSIS:
%   W = addWell(G, rock, W, 'PropertyName', PropertyValue, ...);
%
% PARAMETERS:
%   G       - Grid structure.
%   rock    - Rock data or empty to assume unit permeability
%   W       - Well structure or empty to create new one.
%
%   properties - List of propertyname/propertyvalue pairs.
%                PropertyNames   PropertyValues  Explanation      
%                
%                'Type'      'Dirichlet'       Pressure BCs (faces)
%                            'Neumann'         Flux BCs (faces)
%                            'source'          Sources/sinks (cells)
%                            'bhp'             Pressure constrained well (cells)
%                            'rate'            Rate constrained well (cells)
%
%                'Cells'     <index to cells>
%
%                'Faces'     <index to faces>
%
%                'Values'    <values>          Values / Value corr. to Type
%
%                'Radius'    <value>           Well radius (only for)
%
%                'Direction' 'x','y' or 'z'    Well direction (default 'z')
%
%                
% RETURNS:
%   W       - Updated (or freshly created) well structure, each element
%             of which has the following fields:
%               type           :   well/BC -type ('Dirichlet'/'Neumann'/'source'/'bhp'/'rate')
%               cells / faces  :   Grid cells or faces corresponperforated by this well.
%               values         :   Target control value / values.
%               r              :   Well bore radius.
%               dir            :   Well direction.
%               WI             :   Well index.
%
% SEE ALSO:
%   generateWellSystem.
%
% COMMENTS:

if isempty(W)
   W = struct;
   wellNum = 1;
else
   wellNum = length(W) + 1;
end


props = struct(varargin{:});


if ~isfield(props, 'Type')
   error('No well/BC Type given')
else
   wellType = props.Type
   switch wellType
      case 'Dirichlet'
         if ~isfield(props, 'Faces')
            error('No faces given for Dirichlet BCs')
         else
            W(wellNum).type = 'Dirichlet';
            W(wellNum).faces = props.Faces(:);
            W(wellNum).values = props.Values(:);
         end
      case 'Neumann'
         if ~isfield(props, 'Faces')
            error('No faces given for Neumann BCs')
         else
            W(wellNum).type = 'Neumann';
            W(wellNum).faces = props.Faces(:);
            W(wellNum).values = props.Values(:);
         end
      case 'source'
         if ~isfield(props, 'Cells')
            error('No cells given for sources/sinks')
         else
            W(wellNum).type = 'source';
            W(wellNum).cells = props.Cells(:);
            W(wellNum).values = props.Values(:);
         end
      case 'bhp'
         if ~isfield(props, 'Cells')
            error('No cells given for well')
         else
            W(wellNum).type = 'bhp';
            W(wellNum).cells = props.Cells(:);
            W(wellNum).values = props.Values;
            W(wellNum).r     = props.Radius;
            if ~isfield(props, 'Direction')
               W(wellNum).dir = 'z';
            else
               W(wellNum).dir   = props.Direction;
            end
            W(wellNum).WI    = wellInx(G, rock, W(wellNum));
            W(wellNum).dZ    = getDZ(G, W(wellNum));
         end
      case 'rate'
         if ~isfield(props, 'Cells')
            error('No cells given for well')
         else
            W(wellNum).type = 'rate';
            W(wellNum).cells = props.Cells(:);
            W(wellNum).values = props.Values;
            W(wellNum).r     = props.Radius;
            if ~isfield(props, 'Direction')
               W(wellNum).dir = 'z';
            else
               W(wellNum).dir   = props.Direction;
            end
            W(wellNum).WI    = wellInx(G, rock, W(wellNum));
            W(wellNum).dZ    = getDZ(G, W(wellNum));            
         end
      otherwise
         error(['Unknown well/BC type: ' wellType])
   end
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Private helper functions follow
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


function WI = wellInx(G, rock, well)
[dx, dy, dz] = cellDims(G, well.cells);
if ~isempty(rock)
   k = permDiag(rock, well.cells);
else
   k = ones(length(well.cells), 3);
end

if well.dir == 'x',
   d1 = dy; d2 = dz; l = dx;
   k1 = k(:,2); k2 = k(:,3);
elseif well.dir == 'y',
   d1 = dx; d2 = dz; l = dy;
   k1 = k(:,1); k2 = k(:,3);
elseif well.dir == 'z',
   d1 = dx; d2 = dy; l = dz;
   k1 = k(:,1); k2 = k(:,2);
end

re1 = 0.28 * sqrt((d1.^2).*sqrt(k2 ./ k1) + ...
                  (d2.^2).*sqrt(k1 ./ k2));
re2 = (k2 ./ k1).^(1/4) + (k1 ./ k2).^(1/4);

rw  = well.r;
re  = re1 ./ re2;
ke  = sqrt(k1 .* k2);

WI = 2 * pi * l .* ke ./ log(re / rw);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function dZ = getDZ(G, well);
%dZ = zeros(length(well.cells),1);
centroids = G.cells.centroids(well.cells,:);
[z, ind] = min(centroids(:,3));
[dx, dy, dz] = cellDims(G, well.cells(ind));
bhZ = z - dz/2;
dZ = centroids(:,3) - bhZ;

   
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [dx, dy, dz] = cellDims(G, inx)
n = length(inx);
dx = zeros(n,1); dy = zeros(n,1); dz = zeros(n,1);
faceNo = rldecode((1:G.faces.num)',double(G.faces.numNodes));
cellNo = rldecode((1:G.cells.num)',double(G.cells.numFaces));
for k = 1:n
    % lF = find(cellNo == inx(k));
    % faces = G.cellFaces(lF, 2);
    faces = G.cellFaces(cellNo == inx(k), 1);
    nodes = [];
    for k1 = 1:length(faces)
        % lE = find(faceNo == faces(k1));
        % nodes = unique([nodes; G.faceNodes(lE, 1)]);
        nodes = unique([nodes; G.faceNodes(faceNo == faces(k1), 1)]);
    end
    coords = G.nodes.coords(nodes,:);
    dx(k) = max(coords(:,1)) - min(coords(:,1));
    dy(k) = max(coords(:,2)) - min(coords(:,2));
    dz(k) = max(coords(:,3)) - min(coords(:,3));
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function p = permDiag(rock, inx)
perm = rock.perm(inx,:);
if size(perm, 2) == 1
   p = perm(:, [1, 1, 1]);
elseif size(perm, 2) == 3,
   p = perm;
else
   p = perm(:, [1, 4, 6]);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%}
