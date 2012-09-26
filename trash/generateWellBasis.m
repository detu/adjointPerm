function [basis, rates] = generateWellBasis(G, S, CG, CS, w, weight) %#ok
% generateWellBasis -- Generate basis functions for single well.
%
% SYNOPSIS:
%   [basis, rates] = generateWellBasis(G, S, CG, CS, w, weight)
%
% PARAMETERS:
%   G      - Grid data structure.
%
%   S      - System structure describing the underlying fine grid model,
%            particularly the individual cell flux inner products.
%
%   CG     - Coarse grid structure as defined by function
%            generateCoarseGrid.
%
%   CS     - Coarse system structure as defined by function
%            generateCoarseSystem.
%
%   w      - Individual well from global well structure 'W' as defined by
%            functions addWell and assembleWellSystem.
%
%   weight - Scalar weighting function for generating synthetic basis
%            function source terms.  Sampled in every cell of the
%            underlying fine grid.
%
% RETURNS:
%   basis - sparse matrix of size (numCellFaces-by-numCoarseWellBlocks).
%           Each column representing fluxes for the corresponding
%           block-to-well interface basis.
%
%   rates - sparse matrix of size (numWellCells-by-numCoarseWellBlocks).
%           Each column representing well rates for the corresponding
%           block-to-well interface basis.
%
% SEE ALSO:
%   generateCoarseWellSystem, generateMixedWellBasis

% $Id: generateWellBasis.m 1734 2009-03-16 18:10:56Z bska $

   basis = sparse(S.sizeB(1), length(w.coarseCells));
   rates = sparse(length(w.cells), length(w.coarseCells));
   for k = 1 : length(w.coarseCells),
      %% first reservoir subsystem:
      coarseCell = w.coarseCells(k);
      rIndC = find( CG.cells.subCells(:, coarseCell) );
      rIndB = find( sum(S.C(: ,rIndC), 2) );
      rIndD = find( sum(S.D(rIndB, :), 1) )';

      subBI = S.BI(rIndB, rIndB);
      subC  = S.C(rIndB, rIndC);
      subD  = S.D(rIndB, rIndD);

      subF = zeros(size(rIndB));
      src  = G.cells.volumes(rIndC).*weight(rIndC);
      subG = -src ./ sum(src);
      subH = zeros(size(rIndD));

      %% then well subsystem
      wInd   = find(CG.cells.subCells(w.cells, coarseCell));
      subWBI = w.S.BI(wInd, wInd);
      subWC  = w.S.C(wInd, rIndC);
      subWD  = sparse( length(wInd), 0);

      subWF = zeros(size(wInd));
      subWH = [];

      %% join systems
      BI = blkdiag(subBI, subWBI);
      C  = [subC; subWC];
      D  = blkdiag(subD, subWD);

      f  = [subF; subWF];
      g  = subG;
      h  = [subH; subWH];

      %% solve hybrid system:
      [flux, pres, lam] = schurComplementSymm(BI, C, D, f, g, h);
      Bflux             = [C, D] * [pres; -lam];

      resFlux  = Bflux(1:length(rIndB));
      wellFlux = -flux(length(rIndB)+1:end);

      basis(:,k) = sparse(rIndB, 1, resFlux, S.sizeB(1), 1);
      rates(:,k) = sparse(wInd, 1, wellFlux, length(w.cells), 1);
   end
end
