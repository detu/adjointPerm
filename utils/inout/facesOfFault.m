function [faces] = facesOfFault(field,grdecl,G)
faces =[];
dirstruct = struct('X',0,'Y',1,'Z',2);
if(isfield(grdecl.FAULTS.(field),'multflt'))
    dirnum_vec = grdecl.FAULTS.(field).dir;
    for j=1:length(grdecl.FAULTS.(field).dir)
        dirnum = dirstruct.(dirnum_vec(j));
        %dirnum=zeros(size(dir));
        region = grdecl.FAULTS.(field).cells(j,:);
        if(region(2*dirnum+1) ~= region(2*dirnum+2))
            error(['Wrong fault ind direction ',grdecl.FAULTS.(field).dir(j)]);
        end
        for side=0:0
            region(2*dirnum+1:2*dirnum+2)=region(2*dirnum+1:2*dirnum+2)+side;
            [X,Y,Z] = meshgrid(region(1):region(2),...
                region(3):region(4),...
                region(5):region(6));
            logind =  sub2ind(G.cartDims,X(:),Y(:),Z(:));
            cells =cart2active(G,logind);
            
            if(~isempty(cells))
                cf  = mcolon(double(G.cells.facePos(cells)),double(G.cells.facePos(cells+1)-1));
                if(side == 0)
                    ci = find( G.cellFaces(cf,2) == 2*dirnum + 2 );
                else
                    ci = find( G.cellFaces(cf,2) == 2*dirnum + 1 );
                end
                faces = [faces;G.cellFaces(cf(ci),1)];
                %faces = faces;
                %hfaces = [hfaces;cf(ci)];
            end
        end
    end
    %assert(numel(unique(faces)) == numel(faces));
    faces = unique(faces);
end
