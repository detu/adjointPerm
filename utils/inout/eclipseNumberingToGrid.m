function g_value = eclipseNumberingToGrid(G,e_value,grdecl)
%convert eclipse output to rsol object
%refSol = eclipseOutput2Sol(G,output,grdecl)
g_value = -ones(prod(G.cartDims),1);
ci = find(grdecl.ACTNUM==1);
g_value(ci) = e_value;
g_value = g_value(G.cells.indexMap);
return