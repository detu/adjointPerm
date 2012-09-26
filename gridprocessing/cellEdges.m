function edges = cellEdges(G)
%Extract local-to-global edge numbering for grid cells.
%
% SYNOPSIS:
%   edges = cellEdges(G)
%
% PARAMETERS:
%   G - Grid data structure geometrically discretising a reservoir model.
%
% RETURNS:
%   edges - An m-by-3 array mapping cell numbers to edge numbers.
%           Specifically, if cn(i,1)==j and cn(i,3)==k, then global edge 'k'
%           is one of the edges of cell 'j'.  cn(i,2) describes the position 
%           of the edge in the cell (1 to 12).   
%           In a cornerpoint grid, local edge number of zero indicates that the
%           edge is not part of the original 12 edges. The local numbering
%           corresponds to the edges connecting the following local nodes:
%            1     5
%            2     6
%            3     7
%            4     8
%            1     2
%            3     4
%            5     6
%            7     8
%            1     3
%            2     4
%            5     7
%            6     8
%
% SEE ALSO:
%   grid_structure.

nodes=cellNodes(G);

  edgeNo = rldecode(1:G.cells.num, repmat(12,G.cells.num,1),2).';

  cellToEdges = zeros(numel(edgeNo),2);
  cellToEdges([1:12:end],1) = nodes(find(nodes(:,2)==1),3);
  cellToEdges([1:12:end],2) = nodes(find(nodes(:,2)==5),3);
  cellToEdges([2:12:end],1) = nodes(find(nodes(:,2)==2),3);
  cellToEdges([2:12:end],2) = nodes(find(nodes(:,2)==6),3);
  cellToEdges([3:12:end],1) = nodes(find(nodes(:,2)==3),3);
  cellToEdges([3:12:end],2) = nodes(find(nodes(:,2)==7),3);
  cellToEdges([4:12:end],1) = nodes(find(nodes(:,2)==4),3);
  cellToEdges([4:12:end],2) = nodes(find(nodes(:,2)==8),3);
  cellToEdges([5:12:end],1) = nodes(find(nodes(:,2)==1),3);
  cellToEdges([5:12:end],2) = nodes(find(nodes(:,2)==2),3);
  cellToEdges([6:12:end],1) = nodes(find(nodes(:,2)==3),3);
  cellToEdges([6:12:end],2) = nodes(find(nodes(:,2)==4),3);
  cellToEdges([7:12:end],1) = nodes(find(nodes(:,2)==5),3);
  cellToEdges([7:12:end],2) = nodes(find(nodes(:,2)==6),3);
  cellToEdges([8:12:end],1) = nodes(find(nodes(:,2)==7),3);
  cellToEdges([8:12:end],2) = nodes(find(nodes(:,2)==8),3);
  cellToEdges([9:12:end],1) = nodes(find(nodes(:,2)==1),3);
  cellToEdges([9:12:end],2) = nodes(find(nodes(:,2)==3),3);
  cellToEdges([10:12:end],1) = nodes(find(nodes(:,2)==2),3);
  cellToEdges([10:12:end],2) = nodes(find(nodes(:,2)==4),3);
  cellToEdges([11:12:end],1) = nodes(find(nodes(:,2)==5),3);
  cellToEdges([11:12:end],2) = nodes(find(nodes(:,2)==7),3);
  cellToEdges([12:12:end],1) = nodes(find(nodes(:,2)==6),3);
  cellToEdges([12:12:end],2) = nodes(find(nodes(:,2)==8),3);
  
  [edgesToNodes,i,j] = unique(cellToEdges,'rows');

  edges = zeros(numel(edgeNo),3);
  edges(:,1) = edgeNo;
  edges(:,2) = repmat([1:12]',G.cells.num,1); 
  edges(:,3) = j; 
