function refSol = eclipseOutput2Sol(G,output,grdecl)
%convert eclipse output to rsol object
%refSol = eclipseOutput2Sol(G,output,grdecl)

sat = -ones(prod(G.cartDims),1);
pres = -ones(prod(G.cartDims),1);
ci = find(grdecl.ACTNUM==1);
sat(ci) = output.SWAT.values;
pres(ci) = output.PRESSURE.values;
refSol.s = sat(G.cells.indexMap);
refSol.cellPressure = pres(G.cells.indexMap)*barsa;
return