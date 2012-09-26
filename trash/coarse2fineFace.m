function fineFaces = coarse2fineFace(G, CG, coarseFaces)
% coarse2fineFace -- Convert coarse faces indices to fine face indices 
%                    that are subset of coarse face (e.g. for plotting).
%
% SYNOPSIS:
%   fineFaces = coarse2fineFace(G, CG, coarseFaces)
%
% PARAMETERS:
%   G       - grid_structure data structure.
%
%   CG      - Coarse grid structure
%
%   coarseFaces - Vector of coarse face-indices.
%
% RETURNS:
%   fineFaces - Vector of fine face indices that are subset of coarseFaces.
%
% COMMENTS:
%   not optimized for speed.

% $Id: coarse2fineFace.m 1773 2009-03-19 12:47:11Z bska $

fineFaces = [];
blocks    = double(CG.faces.neighbors(coarseFaces,:));

for i = 1 : size(blocks,1),
   currBlocks = nonzeros(blocks(i,:));
   subCells1  = find(CG.cells.subCells(:,currBlocks(1)));

   %remove faces that are not part of a coarse face
   onCoarseFace = sum(ismember(G.faces.neighbors, subCells1),2) == 1;
   faces = find(onCoarseFace);

   if numel(currBlocks) > 1,  %coarseFace(i) not on boundary
      subCells2 = find(CG.cells.subCells(:,currBlocks(2)));
      fineFInx  = logical(sum(ismember(G.faces.neighbors(onCoarseFace,:),subCells2),2));
   else %boundary face
      fineFInx = G.faces.tag(onCoarseFace) == CG.faces.tag(coarseFaces(i));
   end

   fineFaces = [fineFaces; faces(fineFInx)];
end
