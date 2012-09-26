function [basis, rates] = generateMixedWellBasis(G, S, CG, CS, ...
                                                 w, weight,    ...
                                                 overlapWell,  ...
                                                 overlapBlock) %#ok
% generateMixedWellBasis -- Generate basis functions for single well.
%
% SYNOPSIS:
%   [basis, rates] = generateMixedWellBasis(G, S, CG, CS, w, weight, ...
%                                           overlapW, overlapB)
%
% PARAMETERS:
%   G        - Grid data structure.
%
%   S        - System structure describing the underlying fine grid model,
%              particularly the individual cell flux inner products.
%
%   CG       - Coarse grid structure as defined by function
%              generateCoarseGrid.
%
%   CS       - Coarse system structure as defined by function
%              generateCoarseSystem.
%
%   w        - Individual well from global well structure 'W' as define by
%              functions addWell and assembleWellSystem.
%
%   weight   - Scalar weighting function for generating synthetic
%              basis function source term.  Sampled in every 
%
%   overlapW - Number of fine grid cells by which to extend the support of
%              the well basis function about the well bore.
%
%   overlapB - Number of fine grid cells by which to extend the support of
%              the block-to-well basis function
%
% RETURNS:
%   basis - sparse matrix of size (numCellFaces-by-numCoarseWellBlocks).
%           Each column represents fluxes for the corresponding
%           block-to-well interface basis.
%
%   rates - sparse matrix of size (numWellCells-by-numCoarseWellBlocks).
%           Each column represents well rates for the corresponding
%           block-to-well interface basis.
%
% SEE ALSO:
%   generateCoarseWellSystem, generateMixedBasis.

% $Id: generateMixedWellBasis.m 1820 2009-03-23 12:41:49Z bska $

basis = sparse(S.sizeD(2), length(w.coarseCells));
rates = sparse(length(w.cells), length(w.coarseCells));

% Generate global neighboring connection matrix.
c = (1 : G.cells.num) .';
intF = prod(double(G.faces.neighbors), 2) > 0;
neighborMatrix = sparse([c; double(G.faces.neighbors(intF, 1)); ...
                            double(G.faces.neighbors(intF, 2))],...
                        [c; double(G.faces.neighbors(intF, 2)); ...
                            double(G.faces.neighbors(intF, 1))],...
                        1, G.cells.num, G.cells.num);

% Generate normal sign (sgn) of each cell face.
cellno = rldecode((1 : G.cells.num).', double(G.cells.numFaces));
sgn    = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1) == cellno) - 1;
                     
wellCells = logical( sparse(w.cells, 1, 1, G.cells.num, 1) );
wellOverlap  = []; 
for k = 1: length(w.coarseCells)
   % Define reservoir subsystem component matrices.
   coarseCell   = w.coarseCells(k);
   subCells     = CG.cells.subCells(:, coarseCell);
   
   % Only include wellCells outside block if overlap on well is used. 
   if overlapWell
        wellOverlap  = wellCells;
   end
   blockOverlap = subCells;
   for j = 1 : overlapWell,
      wellOverlap = logical(neighborMatrix * wellOverlap);
   end
   for j = 1 : overlapBlock,
      blockOverlap = logical(neighborMatrix * blockOverlap);
   end
   subCellsOverlap = any([wellOverlap blockOverlap], 2);

   iG0 = subCells;
   iG  = subCellsOverlap;
   iF  = logical(S.C   * iG);
   iH  = logical(S.D.' * iF);

   subB = S.B(iF, iF);
   subC = S.C(iF, iG);
   subD = S.D(iF, iH);

   g0      = zeros([G.cells.num, 1]);
   src     = G.cells.volumes(iG0) .* weight(iG0);
   g0(iG0) = -src ./ sum(src);
   subG    = g0(iG);

   subF = zeros(size(subB, 1), 1);
   subH = zeros(size(subD, 2), 1);

   % Define well subsystem component matrices.
   wInd     = subCellsOverlap(w.cells);
   subWB    = w.S.B(wInd, wInd);
   subWC    = w.S.C(wInd, iG);
   subWD    = ones( nnz(wInd), 1);
   subWDHat = diag(subWD);

   subWF    = zeros(nnz(wInd), 1);
   subWH    = 0;

   % Join reservoir and well subsystem component matrices to form complete
   % system of linear equations.
   B    = blkdiag(subB, subWB   );
   C    = vertcat(subC, subWC   );
   D    = blkdiag(subD, subWD   );
   DHat = blkdiag(subD, subWDHat);

   f    = [subF; subWF];
   g    = subG;
   h    = [subH; subWH];

   orientation       = [sgn(iF); ones(nnz(wInd), 1)];
   n  = length(orientation);   
   Do =  spdiags(orientation, 0, n, n)*DHat;
  
   neumannFaces      = (sum(D) == 1) .';
   neumannFaces(end) = false;   
    
   fv = mixedSymm(B, C, D(:,neumannFaces), f, g, ...
                  h(neumannFaces), Do); %,'Regularize', true);
   resFlux  = fv(1:nnz(iH));
   wellFlux = -fv(nnz(iH)+1:end);

   basis(:,k) = sparse(find(iH), 1, resFlux, S.sizeD(2), 1);
   rates(:,k) = sparse(find(wInd), 1, wellFlux, length(w.cells), 1);
end
