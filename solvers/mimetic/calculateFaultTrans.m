function [G,rock] = calculateFaultTrans(G,grdecl,rock,varargin)
%Compute transmissibilities of faults and add them to grid and optionaly
%change rock properties for multipliers
%
% SYNOPSIS:
%   [G,rock] = calculateFaultTrans(G,grdecl,rock,varargin)
%
% PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   grdecl -  Grdecl struct with eclipse keyword defined
%
%
%   rock    - Rock data structure with valid field 'perm'.  
%OPTIONAL PARAMETERS
%   modify_rock -    modify permeability for cells with transimisibility
%                    multipliers defined that is 
%                    fault multipliers > 1, MULTZ >1 MULTX, MULTY 
%
%   extras_one_grid - put extra informations about faces which are assosited
%                     with fault multipliers > 1, MULTZ >1 MULTX, MULTY 
%
% RETURNS:
%   G - Grid with fault properties 
%   rock - with modified premeabilities
%
% COMMENTS:
%   PLEASE NOTE: mimetic implementaion of 
%                fault multipliers > 1, MULTZ >1 MULTX, MULTY
%                dependent on geological interpretation
%
% SEE ALSO:
%   solveIncompFlowFault

%{
#COPYRIGHT#
%}

opt = struct('modify_rock',false,'extras_on_grid',false)
opt = merge_options(opt, varargin{:});

faces=[];
multflt=[];
if(~isempty(grdecl))
    %%
    fieldnames = fields(grdecl.FAULTS);
    for i = 1:length(fieldnames)
        field = char(fieldnames(i));
        if(isfield(grdecl.FAULTS.(field),'multflt'))
            new_faces = facesOfFault(field,grdecl,G);
            assert(numel(unique(new_faces)) == numel(new_faces))
            if(~isempty(new_faces))
                faces = [faces;new_faces];
                multflt =[multflt;repmat(grdecl.FAULTS.(field).multflt,numel(new_faces),1)];
            end
        end
    end
    %% ad normal faults 
    ci = find(multflt<1.0);          
    G.faultfaces = faces(ci);
    %a bit to much work
    Trans   = computeTrans(G, rock);
    transfaces = 1./accumarray(G.cellFaces(:,1),1./Trans);
    G.faulttrans =  mult2Trans(transfaces,multflt(ci),faces(ci));%transfaces(faces(ci)).*multflt(ci)./(1-multflt(ci));
    %if(opt.extras_on_grid == true)
    G.faultmult = multflt(ci);
    %end
    ff = unique(G.faultfaces);%check that faultsfaces is not dublicated
    %assert(length(ff) == length(G.faulttrans));
    G.faultcellfaces = face2CellFace(G.faultfaces,G);
 
    %% extra information on grid for future applications
    if(opt.extras_on_grid == true)
        % fault multipliers > 1 is treated as multipliers not assosiated with
        % the fault
        ci = find(multflt>1.0);
        %G.faultfaces_cracks = faces(ci);
        %G.faultfaces_cracks_multflt = multflt(ci);
        hf = face2CellFace(faces(ci),G);
        G.cellfacecracks = hf(:,1);
        G.cellfacecracks_mult =  multflt(ci);
        cc = find(hf(:,1) ~= hf(:,2));
        G.cellfacecracks =[G.cellfacecracks;hf(cc,2)];
        G.cellfacecracks_mult = [G.cellfacecracks_mult;multflt(ci(cc))];
    end
    
    %modify rock assuming half of the volume is assoiated with the fault
    if(opt.modify_rock == true)
        ci = find(multflt>1.0);
        cells = G.faces.neighbors(faces(ci));
        if(~isempty(cells))
            hf = face2CellFace(faces(ci),G);
            dir = ceil(G.cellFaces(hf(:,1),2)/2);
            %rock.perm(cells,dir) =rock.perm(cells,dir).*2/(1+1./multflt(ci));
            ind = sub2ind(size(rock.perm),cells,dir);
            rock.perm(ind) = rock.perm(ind).*2./(1+1./multflt(ci));
            
        end
    end

%% finnished with faults     
    
%% add shales to faultfaces MULTZ (MULTX MULTY)
    % same as code in computeTrans
    %vecstruct=[];
    %dirstruct=[];
    
    if(isfield(grdecl,'MULTZ'))
      logcells = find(grdecl.MULTZ < 1);
      if(~isempty(logcells))
          %mult = grdecl.MULTZ(logcells);
          cells = cart2active(G,logcells);
          mult = grdecl.MULTZ(G.cells.indexMap(cells));
          hfaces = mcolon(double(G.cells.facePos(cells)),double(G.cells.facePos(cells+1)-1));
          cellfaces = G.cellFaces(hfaces,:);
          ci = find( cellfaces(:,2) == 6);
          b = cumsum( cellfaces(:,2) == 6);
          hfnum = double(G.cells.facePos(cells+1))-double(G.cells.facePos(cells));
          facepos = cumsum(hfnum);
          numf=b(facepos)-[0;b(facepos(1:end-1))];
          mult = rldecode(mult,numf);
          G.shalefaces = cellfaces(ci,1);
          G.shaletrans = mult2Trans(transfaces,mult,G.shalefaces);
          if(opt.extras_on_grid == true)
              G.shalemult = mult;
          end
          G.shalecellfaces = face2CellFace(G.shalefaces,G);
      end
    end
    
    
    %% add other multipliers not assosiated with faces'
    if(opt.extras_on_grid == true || opt.modify_rock == true)
        vecstruct=[];
        dirstruct=[];
        if(isfield(grdecl,'MULTX'))
            vecstruct = [vecstruct,grdecl.MULTX];
            dirstruct = [dirstruct,2];
            mult = ones(G.cartDims);
            multx = reshape(grdecl.MULTX,G.cartDims);
            mult(2:end,:,:) = multx(1:end-1,:,:);
            vecstruct = [vecstruct,mult(:)];
            dirstruct = [dirstruct,1];%set the direction number
        end
        if(isfield(grdecl,'MULTY'))
            vecstruct = [vecstruct,grdecl.MULTY];
            dirstruct = [dirstruct,4];
            mult = ones(G.cartDims);
            multx = reshape(grdecl.MULTY,G.cartDims);
            mult(:,2:end,:) = multx(:,1:end-1,:);
            vecstruct = [vecstruct,mult(:)];
            dirstruct = [dirstruct,3];
        end
        if(isfield(grdecl,'MULTZ'))
            ci = find(grdecl.MULTZ<1);
            grdecl.MULTZ(ci) = 1;
            vecstruct = [vecstruct,grdecl.MULTZ];
            dirstruct = [dirstruct,6];
            mult = ones(G.cartDims);
            multx = reshape(grdecl.MULTZ,G.cartDims);
            mult(:,:,2:end) = multx(:,:,1:end-1);
            vecstruct = [vecstruct,mult(:)];
            dirstruct = [dirstruct,5];
        end
        
        %modify rock
        if(opt.modify_rock == true)
            if(~isempty(vecstruct))
                for i = 1:length(dirstruct)
                    dim = ceil(dirstruct(i)/2);
                    %
                    rock.perm(:,dim) = rock.perm(:,dim).*2./(1+1./vecstruct(G.cells.indexMap,i));
                end
            end
        end
        if(opt.extras_on_grid == true)
            if(~isempty(vecstruct))
                %% do multiplicators
                for i = 1:length(dirstruct)
                    logcells = find( vecstruct(:,i)>0 & vecstruct(:,i)~= 1);
                    if(~isempty(logcells))
                        mult = vecstruct(logcells,i);
                        dir = dirstruct(i);
                        cells = cart2active(G,logcells);
                        %faces = mcolon(double(G.cells.facePos(cells)),double(G.cells.facePos(cells+1)-1));
                        %cellfaces = G.cellFaces(faces,:);
                        cellfaces = [G.cellFaces,[1:size(G.cellFaces,1)]'];
                        hfaces = mcolon(double(G.cells.facePos(cells)),double(G.cells.facePos(cells+1)-1));
                        cellfaces = cellfaces(hfaces,:);
                        ci = find(cellfaces(:,2) == dir);
                        b = cumsum(cellfaces(:,2) == dir);
                        hfnum = double(G.cells.facePos(cells+1))-double(G.cells.facePos(cells));
                        facepos = cumsum(hfnum);
                        numf=b(facepos)-[0;b(facepos(1:end-1))];
                        mult =rldecode(mult,numf);
                        assert(numel(ci) == numel(mult))
                        G.cellfacecracks = [G.cellfacecracks;cellfaces(ci,3)];
                        G.cellfacecracks_mult = [G.cellfacecracks_mult;mult];
                    end
                    
                end
            end
            
        end
    end
   %% remove fault/shale multiplicators on boundary
   % maybe one shuld remove dublications in faultfaces
   if(isfield(G,'faultfaces'))
    ci = find(G.faces.neighbors(G.faultfaces,2) > 0);
    G.faultfaces = G.faultfaces(ci);
    G.faulttrans = G.faulttrans(ci);
    G.faultmult = G.faultmult(ci);
    G.faultcellfaces =  G.faultcellfaces(ci,:);
   end
   if(isfield(G,'shalefaces'))
    ci = find(G.faces.neighbors(G.shalefaces,2) > 0);
    G.shalefaces = G.shalefaces(ci);
    G.shaletrans = G.shaletrans(ci);
    G.shalecellfaces =  G.shalecellfaces(ci,:);
   end

end
end

function trans = mult2Trans(transfaces,mult,faces)
    trans =  transfaces(faces).*mult./(1-mult);
end
