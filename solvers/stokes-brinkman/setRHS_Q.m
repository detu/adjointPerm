function S = setRHS_Q(S, Pdofs, G, q, type)
% setRHS_Q - Implements values for Q on the right hand side of the
%            Stokes-Brinkman system S
%
% SYNOPSIS:
%   S=setRHS_Q(S, Dofs, G, q, type)
%
% PARAMETERS:
%   S          - The mixed Stokes-Brinkman system
%
%   G          - The grid structure
%
%   Pdofs      - The degrees-of-freedom for the pressure (Dofs.Pdofs)
%
%   q          - Either function handles to calculate the q in each DOF,
%                vectors with values for each DOF for the x-component or a struct.
%
%   type       - The type of q: 'function' for function handles, 'vector' or
%                'ratevector' for vectors of values or 'src' for a struct. String.
% RETURNS:
%   S   - The updated system structure with boundary conditions implemented
%   sol - A structure containing the velocities 'v1', 'v2' and 'v3' and the
%         pressure 'p'
  
  error(nargchk(5, 5, nargin, 'struct'));
  
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
  
  if(~isfield(S,'Q'))
    S.Q = zeros(max(max(Pdofs)),1);
  end
  
  x0 = min(G.nodes.coords(:,1)); 
  y0 = min(G.nodes.coords(:,2));  
 
  dxi = find(G.cellFaces(:,2)==1);  dxj = find(G.cellFaces(:,2)==2);
  dyi = find(G.cellFaces(:,2)==3);  dyj = find(G.cellFaces(:,2)==4);
  dx  = G.faces.centroids(G.cellFaces(dxj,1),1)-...
        G.faces.centroids(G.cellFaces(dxi,1),1); 
  dy  = G.faces.centroids(G.cellFaces(dyj,1),2)-...
        G.faces.centroids(G.cellFaces(dyi,1),2); 
  if dim==3
    z0  = min(G.nodes.coords(:,3));              
    dzi = find(G.cellFaces(:,2)==5);  dzj = find(G.cellFaces(:,2)==6);
    dz  = G.faces.centroids(G.cellFaces(dzj,1),3)-...
          G.faces.centroids(G.cellFaces(dzi,1),3); 
  end

  if dim==3
    [x{1:3}]       = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap); 
    G.cells.ijkMap = double([x{:}]);
  else
    [x{1:2}]       = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap); 
    G.cells.ijkMap = double([x{:}]);
  end

  if(strcmp(type,'function'))
    if(~isa(q,'function_handle'))
      error('q is not a function handle');
    end

    for i=1:G.cells.num,
      % Express x, y and z using xi, eta and zeta
      a = x0+G.cells.ijkMap(i,1)*dx(i)-dx(i)/2;
      b = dx(i)/2;
      c = y0+G.cells.ijkMap(i,2)*dy(i)-dy(i)/2;
      d = dy(i)/2;
      if dim==3
        e = z0+G.cells.ijkMap(i,3)*dz(i)-dz(i)/2;
        f = dz(i)/2;
      end
      
      Pdofs_i = Pdofs(i,:);
      
      if dim==2
        for j=1:length(Pdofs_i),
          refW=L{j};
          integrand1=@(x,y) refW(x,y).*q(a+b*x,c+d*y);
          if strcmp(S.basis, 'BasisFunc2D_CR')
            val1=dblquad(integrand1,-1,1,-1,1)*dx(i)*dy(i);
          else
            val1=dblquad(integrand1,-1,1,-1,1)*dx(i)*dy(i)/4;
          end
          S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+val1;
        end 
      elseif dim==3
        for j=1:length(Pdofs_i),
          refW=L{j};
          integrand1=@(x,y,z) refW(x,y,z).*q(a+b*x,c+d*y,e+f*z);
          if strcmp(S.basis, 'BasisFunc3D_CR')
            val1=triplequad(integrand1,-1,1,-1,1,-1,1)*dx(i)*dy(i)*dz(i);
          else
            val1=triplequad(integrand1,-1,1,-1,1,-1,1)*dx(i)*dy(i)*dz(i)/8;
          end
          S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+val1;
        end 
      end
    end
    
  elseif(strcmp(type,'vector'))
    if(length(q)~=G.cells.num)
      error('q has the wrong size');
    end
    Qind=find(q);
    for kk=1:numel(Qind),
      i=Qind(kk);
      Pdofs_i=Pdofs(i,:);
      for j=1:length(Pdofs_i),
          if strcmp(S.basis, 'BasisFunc2D_CR')
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q(i)*dx(i)*dy(i);
          elseif strcmp(S.basis, 'BasisFunc3D_CR')
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q(i)*dx(i)*dy(i)*dz(i);
          else
            if dim==2
              S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q(i)*dx(i)*dy(i)/4;
            elseif dim==3
              S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q(i)*dx(i)*dy(i)*dz(i)/8;
            end
          end
      end   
    end
    
  elseif(strcmp(type,'ratevector'))
    if(length(q)~=G.cells.num)
      error('q has the wrong size');
    end
    Qind=find(q);
    for kk=1:numel(Qind),
      i=Qind(kk);
      Pdofs_i=Pdofs(i,:);
      for j=1:length(Pdofs_i),
          if strcmp(S.basis, 'BasisFunc2D_CR')
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q.(i);
          else
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q(i)/(2^dim);
          end
      end   
    end
    
  elseif(strcmp(type,'src'))
    for kk=1:numel(q.cell),
      i=q.cell(kk);
      Pdofs_i=Pdofs(i,:);
      for j=1:length(Pdofs_i),
          if strcmp(S.basis, 'BasisFunc2D_CR')
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q.rate(kk);
          else
            S.Q(Pdofs_i(j))=S.Q(Pdofs_i(j))+q.rate(kk)/(2^dim);
          end
      end   
    end

  else
    msg=['Type of RHS_Q can only be ''function'', ''vector'','  ...
         '''ratevector'' or ''src''.'];
    error(msg);
  end
 
