function s = grdecl2Saturation(grdecl,G)
%  convert grdecl structure with water cut to saturation 
%
%  SYNOPSIS
%
%    s = grdecl2Saturation(grdecl,G)
%  
%  PARAMETERS:
%   grdecl - a grdecl structure
%   G   - grid
%
%  RETURN:
%   rock - stucture
%
    s = repmat([1 0 0],G.cells.num,1);    
    if(~isfield(grdecl,'EQLNUM'))
      error('Can not read EQUIL with out EQNUM')
    end
    for i=1:size(grdecl.EQUIL,1)
        logind = find(grdecl.EQLNUM == grdecl.EQUIL(i,7));
        cells = cart2active(G,logind);
        oc = find(G.cells.centroids(cells,3) <= grdecl.EQUIL(i,3) &...
            G.cells.centroids(cells,3) >= grdecl.EQUIL(i,2));
        %s(cells(oc),:) = repmat([0 1 0],length(oc),1);
        s(cells(oc),:) = repmat([0 1 0],length(oc),1);
        oc = find(G.cells.centroids(cells,3) <= grdecl.EQUIL(i,1));
        s(cells(oc),:) = repmat([0 0 1],length(oc),1);
    end