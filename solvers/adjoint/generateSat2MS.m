function [CS] = generateSat2MS(G, S, W, CG, CS, SG, p, fluid)
%
% Create mappings from coarse saturation field to MS - pressure system
% solveMixedWellSystemMS to understand ...

%Lt = fluid.Lt(resSol); % not used here
Lt  = .5*ones( G.cells.num, 1 );
omega = 0.5*ones(G.cells.num,1); 

%[BW, basisW, CW, DW, fW, hW] = unpackWellSystemComponentsMS(W, Lt);
[BW, CW, DW, fW, hW, PsiW, PhiW] = ...
      unpackWellSystemComponentsMS(W, G, p, Lt, omega);   
   
activeCF = CS.activeCellFaces;
activeF  = CS.activeFaces;
numACF   = numel(activeCF);
numAF    = numel(activeF);

cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);

[Bv, Phi] = basisMatrixMixed(G, CG, CS);

BfullBasis = [Bv, -PsiW]; 
fullBasis  = S.BI*BfullBasis; 


h = waitbar(0, 'Computing saturation to MS mappings ...');
for k = 1 : SG.cells.num
   fullBasisDk  = fullBasis' *  spdiags(S.C*SG.cells.subCells(:, k) , 0, size(S.B,1), size(S.B,2));
   %indK         = find( S.C*SG.cells.subCells(:, k) );
   Bk           = fullBasisDk * BfullBasis;
   %Bk           = fullBasisT(:, indK) * BfullBasis( indK, :);
   [i, j, vals] = find( Bk );
   ii{k} = i;
   jj{k} = j;
   vv{k} = vals;
   numVal(k) = length(vals);
   waitbar(k / SG.cells.num, h);
end
close(h)
sat2MS.i = vertcat( ii{:} );
sat2MS.j = vertcat( jj{:} );
i1 = ( 1:length(sat2MS.i) )';
j1 = rldecode( (1:SG.cells.num)', numVal');
sat2MS.vals = sparse(i1, j1, vertcat( vv{:} ));
CS.sat2MS = sat2MS;


% find subCellFaces of coarse cellFaces 
partition    = SG.cells.subCells * (1:SG.cells.num)'; 
cfns         = cellFaceNeighbors(G);
posN         = (cfns(:, 2) > 0);
cfnsSG       = zeros( size(cfns) );
cfnsSG(:, 1)    = partition( cfns(:,1) );
cfnsSG(posN, 2) = partition( cfns(posN, 2) );
subCellFaces    = and( ( cfnsSG(:,1) ~= cfnsSG(:, 2) ), cfnsSG(:, 2) );

CS.subCellFaces = subCellFaces;
ssc             = length( subCellFaces );
%CS.CFFullBasis  = spdiags(subCellFaces, 0, ssc, ssc) * fullBasis;
CS.CFFullBasis  = fullBasis;


function cfns = cellFaceNeighbors(G)
% finds neighbor pairs according to cellface - cells
% [cellface-cells celface-cell-neigbors
cfns     = double( G.faces.neighbors(G.cellFaces(:,1), :) );
cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
sgn      = 2*(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo)-1;
flipRowInx  = ( sgn < 0 );
cfns(flipRowInx, :) = cfns(flipRowInx, [2 1]);
extInd      = find(cfns(:,1)==0);
cfns(extInd, :) = cfns(extInd, [2 1]);
return





