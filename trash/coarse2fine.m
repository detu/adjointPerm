function f = coarse2fine(c, CG)
f = CG.cells.subCells * c;
%f = c(CG.partition)
