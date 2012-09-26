function W = addAdjointWellFields(CG, W, varargin)
%
% SYNOPSIS:
%   W = addAdjointWellFields(CG, W, varargin)
%   W = addAdjointWellFields(CG, W, 'OverlapWell', 0, 'OverlapBlock', 4);                         
%
%
% DESCRIPTION: 
%  Hack for the adjoint code:
%  Add additional fields to the Well-structure.
%  The added fields are coarseCells and optionally CS.overlap and
%  CS.wellOverlap.
%
%  Needed in updateWells and createSingleCellPseudoWells.

opt = struct('OverlapWell', 0, 'OverlapBlock', 0);
opt = merge_options(opt, varargin{:});
overlapW = opt.OverlapWell;
overlapB = opt.OverlapBlock;

numWells = length(W);
for k = 1 : numWells,
   [trash, coarseCells] = find(CG.cells.subCells(W(k).cells,:));
   coarseCells       = unique(coarseCells);
   W(k).coarseCells  = coarseCells;
   
   W(k).CS.overlap     = overlapB;
   W(k).CS.wellOverlap = overlapW;
end
