function T = computeTrans(G, rock, varargin)
%Compute transmissibilities.
%
% SYNOPSIS:
%   T = computeTrans(G, rock)
%
% PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g., perm = perm * darcy() / 1e3 if
%             the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell in 2d/3d and FOUR/SIX columns for a
%             symmetric full tensor permeability.  In the latter case, each
%             cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
% RETURNS:
%   T - half-transmissibilities for each local face of each grid cell
%       in the grid.  The number of half-transmissibilities equal the
%       number of rows in G.cellFaces.
%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas.
%
% SEE ALSO:
%   computeMimeticIP, darcy, permTensor.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-26 15:00:08 +0100 (ma, 26 okt 2009) $
% $Revision: 3076 $

opt = struct('verbose', false,'grdecl',[],'K_system','xyz');
opt = merge_options(opt, varargin{:});

dispif(opt.verbose, 'Computing one-sided transmissibilities...\t');
ticif(opt.verbose);

cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2)';
sgn    = 2*(cellNo == G.faces.neighbors(G.cellFaces(:,1), 1)) -1;
N      = bsxfun(@times, sgn, G.faces.normals(G.cellFaces(:,1),:));
C      = G.faces.centroids(G.cellFaces(:,1), :) - ...
         G.cells.centroids(cellNo, :);
C2     = sum(C.*C, 2);



if(strcmp(opt.K_system,'xyz'))
    [K, i, j] = permTensor(rock, size(G.nodes.coords,2));
    T = sum(C(:,i) .* K(cellNo,:) .* C(:,j), 2) ./ C2 .* ...
        sum(C                     .* N     , 2) ./ C2;
elseif(strcmp(opt.K_system,'loc_xyz'))
    if(size(rock.perm,2) == size(G.cartDims,2))
        dim = ceil(G.cellFaces(:,2)/2);%perm.rock( 
    else
        error('loc_xyz is only valid for diagonal tensor')
    end
    ind = sub2ind(size(rock.perm),cellNo,double(dim));
    T = reshape(rock.perm(ind),numel(ind),1) .* sum(C.* N, 2)./ C2;
else
   error('Unknown coordinate system for permeability'); 
end

nt = find(T<0);
if(length(nt)>0)
    dispif(opt.verbose, ['\nWarning:\n\t', int2str(length(nt)), ...
                         ' negative transmissibilities\n']);
    dispif(opt.verbose, ['\tWe take the absolute values\n']);
    T(nt) = abs(T(nt));
end

% do fault multiplicators which ar in grdecl.FAULTS.('fault_name')
grdecl = opt.grdecl;
dirstruct = struct('X',0,'Y',1,'Z',2);
%trans_mult = ones(size(T));
if(~isempty(opt.grdecl))
    grdecl = opt.grdecl;
    fieldnames = fields(grdecl.FAULTS);
    for i = 1:length(fieldnames)
        field = char(fieldnames(i));
        % field is the fault tag
        if(isfield(grdecl.FAULTS.(field),'multflt'))
            dirnum_vec = grdecl.FAULTS.(field).dir;
            for j=1:length(grdecl.FAULTS.(field).dir)
                dirnum = dirstruct.(dirnum_vec(j));
                %dirnum=zeros(size(dir));
                region = grdecl.FAULTS.(field).cells(j,:);
                if(region(2*dirnum+1) ~= region(2*dirnum+2))
                    error(['Wrong fault ind direction ',grdecl.FAULTS.(field).dir(j)]);
                end
                for side=0:1
                    
                    region(2*dirnum+1:2*dirnum+2)=region(2*dirnum+1:2*dirnum+2)+side;
                    %region(2*dirnum+1:2*dirnum+2)=region(2*dirnum+1:2*dirnum+2)+side;
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
                        cf = cf(ci);
                        T(cf) = T(cf)*grdecl.FAULTS.(field).multflt;
                        %trans_mult(cf) = grdecl.FAULTS.(field).multflt;

                    end
                end
            end
        end
    end    

    %assume MULT?M is not used
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
     vecstruct = [vecstruct,grdecl.MULTZ];
     dirstruct = [dirstruct,6];
     mult = ones(G.cartDims);
     multx = reshape(grdecl.MULTZ,G.cartDims);
     mult(:,:,2:end) = multx(:,:,1:end-1);
     vecstruct = [vecstruct,mult(:)];
     dirstruct = [dirstruct,5];
    end
    
    if(~isempty(vecstruct))
        %% do multiplicators
        for i = 1:length(dirstruct)
            vec = vecstruct(:,i);
            dir = dirstruct(i);
            b = cumsum(G.cellFaces(:,2) == dir);
            cfind = cumsum(double(G.cells.numFaces));%facepos
            numf=b(cfind)-[0;b(cfind(1:end-1))];
            ci = find( G.cellFaces(:,2) == dirstruct(i));
            mult = rldecode(vec(G.cells.indexMap),numf);
            T(ci) = T(ci).*mult;
            %trans_mult(ci) = mult;
        end
    end
    %T = T.*trans_mult;
end
            

tocif(opt.verbose);
