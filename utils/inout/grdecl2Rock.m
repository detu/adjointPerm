function rock = grdecl2Rock(grdecl,grid)
%  convert grdecl structure to rock 
%
%  SYNOPSIS
%
%   rock = grdecl2Rock(grdecl,grid)
%  
%  PARAMETERS:
%   grdecl - a grdecl structure
%   grid   - grid
%
%  RETURN:
%   rock - stucture with permeability and porosity
%         if PERMXY is defiened we assume full tensor permeability 
%
% rock = grdecl2Rock(grdecl,grid)
index = grid.cells.indexMap;
if(isfield(grdecl,'PERMXY'))
    rock.perm = [grdecl.PERMX(index),grdecl.PERMXY(index),grdecl.PERMZX(index),...
        grdecl.PERMY(index),grdecl.PERMYZ(index),grdecl.PERMZ(index)];
    rock.poro = grdecl.PORO(index);
else
    rock.perm = [grdecl.PERMX(index),grdecl.PERMX(index),grdecl.PERMZ(index)];
    rock.poro = grdecl.PORO(index);
end
rock.perm = rock.perm*milli*darcy;
    
end

