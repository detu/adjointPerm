function SG = generateNonUniformCoarseGridv2(G, solution, rock, S, lL, ...
                                           lU, printfigs, type)
% Hacked from Veras version .....
%
% generateNonUniformCoarseGrid -- generates non-uniform coarse grid
% according to Jørg Aarnes strategy, described in "Coarsening of
% three-dimensional structured and unstructured grids for subsurface
% flow" [mer ref.?]. Follows the four steps:
% Step 1: Initiate the coarse grid according to magnitude of the
% velocity.
% Step 2: Merge too small blocks with neighboring blocks.
% Step 3: Refine too large blocks into smaller blocks.
% Step 4: Merge too small blocks again.
% Coarsening parameters:
%  lL - lower bound on volume of blocks
%  lU - upper bound on total amount of flow in blocks.
%
% SYNOPSIS:
%   CG = generateNonUniformCoarseGrid(G,solution, rock, S, BC, lL, lU)
%
% PARAMETERS:
%  G         - grid structure
%  solution  - stucture of the solution of the linears system S
%  rock      - structure with grid data
%  S         - structure of the linear system
%  BC        - structure of the boundary conditions
%  lL        - coarsening parameter, giving the lower bound on volume of blocks
%  lU        - coarsening parameter, giving the upper bound on total
%              amount of flow in blocks.
%
% RETURNS:
%  CG        - structure of the coarse grid. Additional data is added by
%              call of function simulationGridData
%
% REFERENCE:
%   See also: generateCartesianCoarseGrid, computeCellFlux, initSimGrid,
%   postprocessGrid, refineGrid, simulationGridData,

%
% TODO:
%  Remove counting of cells,
%  option of plotting

% $Id: generateNonUniformCoarseGrid.m 788 2008-09-29 08:56:59Z vlh $
addpath([SAMSIMROOTDIR, filesep, 'examples', filesep, 'simpleNUC']) 

cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);

q = S.C'*solution.cellFlux;

sourceCells = q.*(q > 0); sinkCells = q.*(q < 0);

%% Compute a representative flux value for the cell center.
cellFlux = computeCellFlux(G, solution.faceFlux./G.faces.areas, ...
                           sourceCells, sinkCells);
if printfigs,
   figure
   plotCellData(G, log10(cellFlux));
   title('Cell Flux')
   axis tight
   colorbar
end

%% Grid coarsening
%% Step 1:
CG = initSimGrid(G, cellFlux); CG.num
if printfigs,
   figure
   subplot(2,2,1)
   PlotDgrid(CG, G);
   axis tight
end

%CG=simulationGridData(CG,G, rock, S.RHS.g, inflow, outflow);
%Fluxes=simFlux(CG,G,solution.faceFlux,inflow,outflow);
%sum(sum(spones(CG.coarseGrid)))

%% Step 2:
CG = postprocessGrid(CG, G, lL, log10(cellFlux)); CG.num
if printfigs,
   subplot(2,2,2)
   PlotDgrid(CG, G);
   axis tight
end
%sum(sum(spones(CG.coarseGrid)))

%% Step 3:
CG = refineGrid(CG, G, lU, log10(cellFlux), type); CG.num
if printfigs,
   subplot(2,2,3)
   PlotDgrid(CG,G);
   axis tight
end
%sum(sum(spones(CG.coarseGrid)))

%% Step 4:
CG = postprocessGrid(CG, G, lL, log10(cellFlux)); CG.num
if printfigs,
   subplot(2,2,4)
   PlotDgrid(CG,G);
   axis tight
end

SG.cells.num  = CG.num;
SG.cells.subCells = sparse(G.cells.num, SG.cells.num);
for k = 1 : SG.cells.num
    [i,j,v] = find( CG.coarseGrid(:,k) );
    SG.cells.subCells(v, k) = 1;
end
    
%sum(sum(spones(CG.coarseGrid)))
return
CG = simulationGridData(CG, G, rock, S.RHS.g_src, ...
                        [], []);

Fluxes = simFlux(CG, G, solution.faceFlux, ...
                 [], []);

sum(sum(spones(CG.coarseGrid)))
