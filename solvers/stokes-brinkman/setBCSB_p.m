function S = setBCSB_p(S,BC,Dofs,G, varargin) 
% setBCSB_p - Implements the pressure boundary conditions into the system S
%
% SYNOPSIS:
%   S = setBCSB_p(S,BC,Dofs,G) 
%
% PARAMETERS:
%   S          - The mixed Stokes-Brinkman system
%
%   BC         - Boundary condition structure from a prior call to 'addBCSB' or
%                an empty structure, which will automatically give no-flow conditions.
%
%   Dofs       - The degrees-of-freedom structure
%
%   G          - The grid structure
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - dx -- A structure with fields 'dx', 'dy', (and 'dz')
%                 containing a vector with values for each cell.
%
% RETURNS:
%   S   - The updated system structure with boundary conditions implemented
%

  error(nargchk(4, 6, nargin, 'struct'));
  opt   = struct('dx', []);
  opt   = merge_options(opt, varargin{:});
  dxvec = opt.dx;
  
  dim = numel(G.cartDims);
  
  % Read the basis functions
  if strcmp(S.basis, 'TH')
    if dim==2
      load 'BasisFunc2D';
    elseif dim==3
      load 'BasisFunc3D';
    end
  else
    load(S.basis);
  end
  
  % Choose the pressure boundary conditions
  isP    = BC.nodetype == 2; %strcmp('2',BC.nodetype);
  BCdofs = BC.node(isP);
  BCval  = BC.nodevalue(isP);
  BCtag  = BC.nodetag(isP);
  
  % Separate the upper, lower, left and right boundaries
  ile = find(BCtag==1); BCdofs_le = unique(BCdofs(ile)); BCval_le = BCval(ile);
  ir  = find(BCtag==2); BCdofs_r  = unique(BCdofs(ir));  BCval_r  = BCval(ir);
  il  = find(BCtag==3); BCdofs_l  = unique(BCdofs(il));  BCval_l  = BCval(il);
  iu  = find(BCtag==4); BCdofs_u  = unique(BCdofs(iu));  BCval_u  = BCval(iu);
  ib  = find(BCtag==5); BCdofs_b  = unique(BCdofs(ib));  BCval_b  = BCval(ib);
  it  = find(BCtag==6); BCdofs_t  = unique(BCdofs(it));  BCval_t  = BCval(it);

  if isempty(dxvec)
    dxi = find(G.cellFaces(:,2)==1);  dxj = find(G.cellFaces(:,2)==2);
    dyi = find(G.cellFaces(:,2)==3);  dyj = find(G.cellFaces(:,2)==4);
    dx  = G.faces.centroids(G.cellFaces(dxj,1),1)-...
          G.faces.centroids(G.cellFaces(dxi,1),1);
    dy  = G.faces.centroids(G.cellFaces(dyj,1),2)-...
          G.faces.centroids(G.cellFaces(dyi,1),2);
    if dim==3
      dzi = find(G.cellFaces(:,2)==5);  dzj = find(G.cellFaces(:,2)==6);
      dz  = G.faces.centroids(G.cellFaces(dzj,1),3)- ...
            G.faces.centroids(G.cellFaces(dzi,1),3); 
    end
  else
    dx = dxvec.dx;
    dy = dxvec.dy;
    if dim==3
      dz = dxvec.dz;
    end
  end

  % Modify S.F1 and S.F2 

  for i=1:G.cells.num,
     if dim==2
      % Integrate over the left boundary
      if(~isempty(BCdofs_le))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_le);
        for j = 1:length(IV),
          refW             = W_1y{IV(j)};
          S.F1(Vd_kant(j)) = S.F1(Vd_kant(j))+...
              refW*BCval_le(find(BCdofs_le==Vd_kant(j)))*dy(i)/2;
             refW*BCval_le(find(BCdofs_le==Vd_kant(j)))*dy(i)/2;
        end
      end

      % Integrate over the right boundary
      if(~isempty(BCdofs_r))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_r);
        for j = 1:length(IV),
          refW             = W1y{IV(j)};
          S.F1(Vd_kant(j)) = S.F1(Vd_kant(j))-...
              refW*BCval_r(find(BCdofs_r==Vd_kant(j)))*dy(i)/2;
        end,
      end,
      
      % Integrate over the lower boundary
      if(~isempty(BCdofs_l))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_l);
        for j = 1:length(IV),
          refW             = Wx_1{IV(j)};
          S.F2(Vd_kant(j)) = S.F2(Vd_kant(j))+...
              refW*BCval_l(find(BCdofs_l==Vd_kant(j)))*dx(i)/2;
        end
      end  

      % Integrate over the upper boundary
      if(~isempty(BCdofs_u))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_u);
        for j = 1:length(IV),
          refW             = Wx1{IV(j)};
          S.F2(Vd_kant(j)) = S.F2(Vd_kant(j))-...
              refW*BCval_u(find(BCdofs_u==Vd_kant(j)))*dx(i)/2;
        end
      end
 
    elseif dim==3
      
      % Integrate over the left boundary
      if(~isempty(BCdofs_le))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_le);
        for j = 1:length(IV),
          refW             = W_1yz{IV(j)};
          S.F1(Vd_kant(j)) = S.F1(Vd_kant(j))+...
                             refW*BCval_le(find(BCdofs_le==Vd_kant(j)))*dy(i)*dz(i)/4;
        end
      end
      
      % Integrate over the right boundary
      if(~isempty(BCdofs_r))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_r);
        for j = 1:length(IV),
          refW             = W1yz{IV(j)};
          S.F1(Vd_kant(j)) = S.F1(Vd_kant(j))-...
                             refW*BCval_r(find(BCdofs_r==Vd_kant(j)))*dy(i)*dz(i)/4;
        end
      end
      
      % Integrate over the lower boundary
      if(~isempty(BCdofs_l))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_l);
        for j = 1:length(IV),
          refW             = Wx_1z{IV(j)};
          S.F2(Vd_kant(j)) = S.F2(Vd_kant(j))+...
                             refW*BCval_l(find(BCdofs_l==Vd_kant(j)))*dx(i)*dz(i)/4;
        end
      end  
      
      % Integrate over the upper boundary
      if(~isempty(BCdofs_u))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_u);
        for j = 1:length(IV),
          refW             = Wx1z{IV(j)};
          S.F2(Vd_kant(j)) = S.F2(Vd_kant(j))-...
                             refW*BCval_u(find(BCdofs_u==Vd_kant(j)))*dx(i)*dz(i)/4;
        end
      end
      
      % Integrate over the bottom boundary
      if(~isempty(BCdofs_b))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_b);
        for j = 1:length(IV),
          refW             = Wxy_1{IV(j)};
          S.F3(Vd_kant(j)) = S.F3(Vd_kant(j))+...
                             refW*BCval_b(find(BCdofs_b==Vd_kant(j)))*dy(i)*dx(i)/4;
                         
        end
      end
      
      % Integrate over the top boundary
      if(~isempty(BCdofs_t))
        [Vd_kant,IV] = intersect(Dofs.Vdofs(i,:),BCdofs_t);
        for j = 1:length(IV),
          refW             = Wxy1{IV(j)};
          S.F3(Vd_kant(j)) = S.F3(Vd_kant(j))-...
                             refW*BCval_t(find(BCdofs_t==Vd_kant(j)))*dy(i)*dx(i)/4;              
        end 
      end 
     end
   end
  
