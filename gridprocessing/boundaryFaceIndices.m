function ix = boundaryFaceIndices(G, direction, i1, i2, i3)
%Retrieve face indices belonging to subset of global outer faces.
%
% SYNOPSIS:
%   ix = boundaryFaceIndices(G, side, i1, i2, i3)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   side  - Global side from which to extract face indices.
%           Must be one of {'LEFT', 'RIGHT', 'FRONT', 'BACK', ...
%                           'BOTTOM', 'TOP'}
%
%   i1,i2 - Index ranges for local (in-plane) axes one and two,
%           respectively.
%
%   i3    - Index range for global axis perpendicular to 'side'.  The
%           primary purpose of this parameter is to exclude faces *within*
%           a reservoir from being added to the return value 'ix'.  In
%           other words, to have the return value 'ix' only reflect true
%           boundary faces.  This will typically occur in faulted
%           reservoirs where a given face may considered external by virtue
%           of being connected to only a single reservoir cell.
%
% RETURNS:
%   ix    - Required face indices.
%
% NOTE:
%   This function is mainly intended for internal use by functions fluxside
%   and pside.  Its calling interface may change more frequently than those
%   of the *side functions.
%
% SEE ALSO:
%   fluxside, pside.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

st = dbstack(1);
try
   caller = st(1).name;
catch  %#ok
   caller = 'BASE:BoundaryConditions';
end

if isempty(i1) || isempty(i2),
   error([caller, ':IndexSpec:Empty'], ...
         'Empty local index range for boundary condition.');
end

i1 = reshape(i1, [], 1);
i2 = reshape(i2, [], 1);
i3 = reshape(i3, [], 1);

%% Extract all faces of cells within the given subset
[cells, ft, isOutF] = bdryCellsSubset(G, direction, i1, i2, i3, caller);

fIX   = cumsum([0; double(G.cells.numFaces)]);
hfIX  = mcolon(fIX(cells) + 1, fIX(cells + 1));
faces = G.cellFaces(hfIX, 1);
tags  = G.cellFaces(hfIX, 2);

%% Extract only those faces which have the required tag
ix = faces(isOutF(faces) & tags == ft);


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function [cells, faceTag, isOutF] = ...
      bdryCellsSubset(G, direction, i1, i2, i3, caller)
%% Determine which indices and face tags to look for
switch upper(direction),
   case 'LEFT',
      % I == min(I)
      [d1, d2, d3, faceTag] = deal(2, 3, 1, 1);
   case 'RIGHT',
      % I == max(I)
      [d1, d2, d3, faceTag] = deal(2, 3, 1, 2);
   case 'BACK',
      % J == min(J)
      [d1, d2, d3, faceTag] = deal(1, 3, 2, 3);
   case 'FRONT',
      % J == max(J)
      [d1, d2, d3, faceTag] = deal(1, 3, 2, 4);
   case 'TOP',
      % K == min(K)
      [d1, d2, d3, faceTag] = deal(1, 2, 3, 5);
   case 'BOTTOM',
      % K == max(K)
      [d1, d2, d3, faceTag] = deal(1, 2, 3, 6);
   otherwise
      error([caller, ':Side:Unknown'],                  ...
            'Boundary side ''%s'' not supported in %s', ...
            direction, caller);
end

%% Determine uniqe outer cells (i.e., cells on boundary)
isOutF = any(double(G.faces.neighbors) == 0, 2);
cells  = find(accumarray(sum(G.faces.neighbors(isOutF,:), 2), 1) > 0);

%% Determine logical indices of these cells
% Assume we will only be called for logically Cartesian grids for which the
% fields 'G.cartDims' and 'G.cells.indexMap' are present.
%
dims = reshape(G.cartDims, 1, []);
if size(G.nodes.coords,2) == 2,
   dims = [dims(1:2), 1];

   if faceTag > 4,
      error([caller, ':Side:Unsupported'], ...
            ['Boundary side ''%s'' is not defined for two-dimensional ', ...
             'grids in ''%s''.'], direction, caller);
   end

   if numel(i2) > 1,
      error([caller, ':I2:ERANGE'], ...
            ['Two-dimensional boundary faces are incompatible with a\n', ...
             'two-dimensional grid model in ''%s''.\n\n'               , ...
             'Specifically, ''I2'' must contain only a single number.'], ...
            caller);
   end

   if numel(i3) > 0 && any(i3 > 1),
      error([caller, ':I3:ERANGE'], ...
            ['A non-zero cell depth is incompatible with a\n', ...
             'two-dimensional grid model in ''%s''.\n\n'     , ...
             'Specifically, ''I3'' must be empty.'], caller);
   end
end
[cIJK{1:3}] = ind2sub(dims, double(G.cells.indexMap(cells)));

%% Determine whether or not a given cell is within the required subset
%
if any(i1 < 1 | dims(d1) < i1),
   error([caller, ':I1:ERANGE'], ...
          'Cell range ''I1'' outside model in ''%s''.', caller);
end
if any(i2 < 1 | dims(d2) < i2),
   error([caller, ':I2:ERANGE'], ...
          'Cell range ''I2'' outside model in ''%s''.', caller);
end
if numel(i3) > 0 && any(i3 < 1 | dims(d3) < i3),
   error([caller, ':I3:ERANGE'], ...
          'Cell range ''I3'' outside model in ''%s''.', caller);
end

[I{1:2}] = deal(false([G.cells.num, 1]));
I{1}(i1) = true;
I{2}(i2) = true;
inSubSet = I{1}(cIJK{d1}) & I{2}(cIJK{d2});

if ~isempty(i3),
   I{3}     = false([G.cells.num, 1]);   I{3}(i3) = true;
   inSubSet = inSubSet & I{3}(cIJK{d3});
end

%% Extract required cell subset
cells = cells(inSubSet);
