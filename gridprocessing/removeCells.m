function [G, cellmap, facemap, nodemap] = removeCells(G, cells)
%Remove cells from grid and renumber cells, faces and nodes.
%
% SYNOPSIS:
%   G = removeCells(G, cells);
%   [G, cellmap, facemap, nodemap] = removeCells(G, cells)
%
% PARAMETERS:
%   G          - Valid grid definition
%
%   cells      - list of cell numbers to be removed.
%
%
% RETURNS:
%   G          - Updated grid definition where cells have been removed.
%                In addition, any unreferenced faces or nodes are
%                subsequently removed.
%
% EXAMPLE:
%
%      G = cartGrid([3,5,7]);
%      G = removeCells(G, (1:2:G.cells.num));
%      plotGrid(G);view(-35,20);camlight
%
%
%
% NOTE:
%
%   The process of removing cells is irreversible.
%
% SEE ALSO:
%   readGRDECL, buildMatchingGrid, deactivateZeroPoro.

%{
#COPYRIGHT#
%}

% $Id: removeCells.m 2115 2009-04-29 14:02:36Z bska $

%
% Mark resulting outer boundary with optional tag.
%
%
% If we make tag,
%  m1 = ismember(G.faces.neighbors(:,1), cells);
%  m2 = ismember(G.faces.neighbors(:,2), cells);
%  t  = find((m1 & ~m2) | (~m1 & m2));
%  G.faces.tag(t)=99;


  % New numbering of cells
  ind        = false(G.cells.num,1);
  ind(cells) = true;
  cellmap    = mapExcluding(ind);

  % remove and renumber cells in cellFaces
  G.cellFaces(rldecode(ind, double(G.cells.numFaces)), :) = [];

  % Alter cell numbering in G.faces.neighbors
  n = G.faces.neighbors;
  G.faces.neighbors(n(:,1)>0,1) = cellmap(n(n(:,1)>0,1));
  G.faces.neighbors(n(:,2)>0,2) = cellmap(n(n(:,2)>0,2));
  clear n

  % Alter cells
  G.cells.num             = G.cells.num-numel(cells);
  G.cells.numFaces(cells) = [];
  G.cells.indexMap(cells) = [];


  %new numbering of faces.
  ind     = all(G.faces.neighbors(:,1:2)==0,2);
  facemap = mapExcluding(ind);

  % remove and renumber faces in faceNodes
  G.faceNodes(rldecode(ind, double(G.faces.numNodes))) = [];

  % remove and renumber faces in cellFaces
  G.cellFaces(:,1) = facemap(G.cellFaces(:,1));


  if any(G.cellFaces(:,1)==0),
      error('In removeCells: Too many faces removed!');
  end


  G.faces.neighbors(ind,:) = [];
  G.faces.numNodes (ind,:) = [];
  G.faces.tag      (ind,:) = [];
  G.faces.num              = G.faces.num - sum(ind);



  % Construct node map:
  ind = true(G.nodes.num, 1);
  ind(G.faceNodes) = false;
  nodemap = mapExcluding(ind);

  % Remove nodes
  G.nodes.coords(ind,:) = [];
  G.nodes.num           = G.nodes.num - sum(ind);
  G.faceNodes           = uint32(nodemap(G.faceNodes));

  if any(G.faceNodes==0),
      error('In removeCells: Too many nodes removed!');
  end

%  G.cells.actnum(cells)=false;

  function m = mapExcluding(indices)
      n            = numel(indices);
      ind          = ones(n,1);
      ind(indices) = 0;
      m            = cumsum(ind);
      m(indices)   = 0;
