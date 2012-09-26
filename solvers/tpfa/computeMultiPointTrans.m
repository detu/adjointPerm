function T = computeMultiPointTrans(g, rock)
%Compute multi-point transmissibilities.
%
% SYNOPSIS:
%   T = computeMultiPointTrans(G, rock)
%
% PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and FOUR/SIX columns for
%             a symmetric full tensor permeability.  In the latter case,
%             each cell gets the permeability tensor
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

%{
#COPYRIGHT#
%}

% $Date: 2009-10-14 12:10:31 +0200 (on, 14 okt 2009) $
% $Revision: 3002 $

cellNo   =  rldecode(1:g.cells.num, double(g.cells.numFaces), 2) .'; 
%sgn      =  2*cellNo == g.faces.neighbors(g.cellFaces(:,1),1) - 1;
cells    =  rldecode(g.faces.neighbors, double(g.faces.numNodes));
nodes    =  repmat(g.faceNodes, [2,1]);
faces    =  repmat(rldecode(1:g.faces.num, double(g.faces.numNodes),2)', [2,1]);
subfaces =  repmat((1:size(g.faceNodes,1))', [2,1]);
i        =  cells~=0;
w        =  double(sortrows([cells(i), nodes(i), faces(i), subfaces(i)]));
clearvars cells nodes faces subfaces i

% Make a sparse matrix: use b as block size.
[a, blocksz] = rlencode(w(:,[1,2])); 
dims = size(g.nodes.coords, 2);
assert(all(blocksz==dims));

[i,j] = ndgrid(1:size(w,1), 1:dims);
j     = bsxfun(@plus, j, rldecode(cumsum(blocksz)-blocksz(1), blocksz));

R     = g.faces.centroids(w(:,3),:) - g.cells.centroids(w(:,1),:);
R     = sparse(i,j,R);


sgn   = 2*(w(:,1) == g.faces.neighbors(w(:,3),1)) -1;
N     = g.faces.normals  (w(:,3),:).*sgn(:,ones(1, dims))./ ...
        double(g.faces.numNodes(w(:,3), ones(1,dims)));
N     = sparse(i,j,N);

k     = permTensor(rock, dims);
K     = sparse(i,j,reshape(k(a(:,1),:)', dims, [])');
clearvars i j k a blocksz sgn


iB    = N*K/R;
B     = iB\speye(size(iB));

nhf   = size(w, 1);
e     = (1:nhf)';
s     = 2*(w(:,1)==g.faces.neighbors(w(:,3),1))-1;
D     = sparse(e, w(:,4), 1);
Dm    = sparse(e, w(:,4), s);
C     = sparse(e, w(:,1), 1);
clearvars N R K

DmBDm = Dm'*B*Dm;

tic
[a,b,c,d]=dmperm(DmBDm);
if all(diff(c)==1), % diagonal innerproduct
   iDmBDm = inv(DmBDm);
else
   sz = sum(accumarray(g.faceNodes, 1).^2);
   I = zeros(sz, 1);
   J = I;
   V = I;
   pos = 0;
   for k=1:numel(c)-1,
      i     = a( c(k) : c(k+1)-1 );
      j     = b( d(k) : d(k+1)-1 );
      ip    = inv(full(DmBDm(i,j)));
      sz    = numel(i)^2;

      I(pos+1:pos+sz)  = repmat(i, [numel(j), 1])';
      J(pos+1:pos+sz)  = repmat(j, [numel(i), 1]);
      V(pos+1:pos+sz)  = ip;

      pos   = pos + sz;
   end
   clear DmBDm
   iDmBDm = sparse(I,J,V);
end
toc
tic
%}

%map = sparse(g.cells.num*g.faces.num, 1);
e      = (1:size(g.cellFaces,1))';

% Construct map from cell and subface to half-face
map    = sparse(cellNo*g.faces.num+ double(g.cellFaces(:,1)), 1, e);

% c1 adds up sub-half-face contributions for each half-face
c1     = sparse(1:size(w, 1), map(w(:,1)*g.faces.num + w(:,3)), 1);

% d1 adds up sub-face contributions for each face.
d1     = sparse(w(:,4), w(:,3), 1);

% Compute multi-point transmissibilities for all faces in terms
% of cell pressures and outer boundary pressures.  
b = full(sum(D, 1)) == 1;
%T = c1'*Dm*inv(Dm'*B*Dm)*Dm'*[C, -D(:,b)*d1(b,:)];
T = c1'*Dm*iDmBDm*Dm'*[C, -D(:,b)*d1(b,:)];

% c1'*D*d1 har samme struktur som vanlig D.
% T er feil størrelse å returnere dersom gravitasjon skal håndteres
% skikkelig.  Gravitasjonsleddet kommer inn som c1'*Dm*iDmBDm*Dm'*f.
% Siden gravitasjonsbidraget for all subfaces er likt kan de sikkert
% skrives om til c1'*Dm*iDmBDm*Dm'*F*g der F*G=f.

%T(:,all(T==0, 1))=[];
