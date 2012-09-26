function h = plotInterpCellData(G, data)
%Plot exterior grid faces, coloured by interpolated data, to current axes.
%
% SYNOPSIS:
%   h = plotCellData(G, data);
%   h = plotCellData(G, data, subDomain);
%   h = plotCellData(G, data, subDomain, 'pn1', pv1, ...);
%
% PARAMETERS:
%   G         - Grid data structure.
%
%   data      - Scalar cell data with which to colour the grid.  One scalar
%               value for each cell in the grid (i.e. always a vector of
%               length G.cells.num, even when NUMEL(subDomain) <
%               G.cells.num).
%
%   subDomain - Vector of cell-indices defining sub-grid. OPTIONAL.
%               If unspecified, or empty, the entire grid is coloured (i.e.
%               ISEMPTY(subDomain) -> subDomain = (1 : G.cells.num).').
%
%   'pn'/pv   - List of property names/property values. OPTIONAL.
%               Default value = {'FaceColor', 'flat', 'EdgeColor', 'none'}.
%               All other PATCH property specifications are allowed.
%
% RETURNS:
%   h - Structure of handles to topological unique patch-objects,
%       i.e. h(1) refers to all triangles, h(2) to all quads &c.
%       The structure has the following fields:
%         - patch   :  handle to current patch object
%         - polyN   :  current polygon-type
%
% SEE ALSO:
%   plotCellData, plotGrid, patch.

%{
#COPYRIGHT#
%}

% $Id: plotInterpCellData.m 1952 2009-03-31 10:30:47Z bska $

[f, v, fvc] = getFaceVertexData(G, data);

h = patch('Faces', f, 'Vertices', v, 'FaceVertexCData', fvc, ...
          'FaceColor', 'interp', 'EdgeColor', 'interp', ...
          'FaceLighting', 'phong', 'EdgeLighting', 'phong');

set(gca, 'ZDir', 'reverse')

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [f, v, fvc] = getFaceVertexData(G, data)
% Compute vertices
[v, trash, renum] = unique(G.nodes.coords(:, [3, 2, 1]), 'rows');
v = v(:, [3, 2, 1]);

% Compute faces
fn = renum(G.faceNodes(:,1));

m = G.faces.num;
n = double(max(G.faces.numNodes));
f = nan([n, m]);
o = reshape((0 : m-1) .* n, [], 1);
f(mcolon(o + 1, o + double(G.faces.numNodes))) = fn;
f = f .';  % faces-by-vertices

% Compute FaceVertexCData
% One value per vertex for INTERP colouring/shading (indexed colour system)
eIX = cumsum([0; double(G.faces.numNodes)]);
cf  = rldecode((1 : G.cells.num) .', double(G.cells.numFaces));
nn  = accumarray(cf, G.faces.numNodes(G.cellFaces(:,1)));
cn  = fn(mcolon(eIX(G.cellFaces(:,1)) + 1, ...
                eIX(G.cellFaces(:,1) + 1)));
fvc = accumarray(cn, data(rldecode((1 : G.cells.num).', double(nn)))) ./ ...
      accumarray(cn, 1);
