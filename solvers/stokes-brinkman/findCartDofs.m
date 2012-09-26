function Dofs = findCartDofs(G)
% findCartDofs - Construct 2/3-D degrees-of-freedom structure for Taylor-Hood elements
%
% SYNOPSIS:
%   Dofs = findCartDofs(G)
%
% DESCRIPTION:
%
%  In 2 dimensions:
%  
%   The DOFs are numbered in the following way for a 2 x 2 system:
%
%      Velocity:           Pressure:
%
%   7--24--8--25--9       7----8----9
%   |      |      |       |    |    |
%  17  12 18  13  19      |    |    |
%   |      |      |       |    |    |
%   4--22--5--23--6       4----5----6
%   |      |      |       |    |    |
%  14  10 15  11  16      |    |    |
%   |      |      |       |    |    |
%   1--20--2--21--3       1----2----3
%
%   The local velocity DOFs for each cell are numbered in the following way: 
%
%      Velocity:            Pressure:
%
%       7--8--9              3-----4
%       |  |  |              |  |  |
%       4--5--6              |-----|
%       |  |  |              |  |  |
%       1--2--3              1-----2
%
%  In 3 dimensions:
%  
%   The DOFs are numbered in the following way for a 2 x 2 x 1 system:
%
%  7|--56--8--57--|9         49|--33-50--34--|51        16|--62-17--63--|18
%   |      |      |            |      |      |            |      |      |
% 67|  37 68  38  |69        26|  21 27  22  |28        73|  41 74  42  |75 
%   |      |      |   for      |      |      |   for      |      |      |   for 
%  4|--54--5--55--|6  the    46|--31-47--32--|48 the    13|--60-14--61--|15 the
%   |      |      |   lower    |      |      |   middle   |      |      |   upper 
% 64|  35 65  36  |66 layer, 23|  19 24  20  |25 layer, 70|  39 71  40  |72 layer.
%   |      |      |            |      |      |            |      |      |
%  1|--52--2--53--|3         43|--29-44--30--|45        10|--58-11--59--|12
%
%   The local DOFs for each cell are numbered in the following way: 
%
%   7--8--9  16--17--18  25--26--27
%   |  |  |   |  |   |    |  |   |
%   4--5--6  13--14--15  22--23--24 
%   |  |  |   |  |   |    |  |   |
%   1--2--3  10--11--12  19--20--21
%
% PARAMETERS:
%   G - Grid structure 
%  
% RETURNS:
%   Dofs - degrees-of-freedom structure for Taylor-Hood elements
%          having the following fields: 
%            - Pdofs           :DOFs for the pressure (in the vertices)
%            - numPdofs        :the total number of Pdofs (Nx+1)*(Ny+1)
%                               in 2-D and (Nx+1)*(Ny+1)*(Nz+1) in 3-D
%            - Vdofs           :DOFs for the velocity (in the vertices,
%                               cell centers and edge centers)
%            - numVdofs        :the total number of Vdofs
%                               (2*Ny+1)*(2*Nx+1) in 2-D and
%                               (2*Ny+1)*(2*Nx+1)*(2*Nz+1) in 3-D 
%            - HalfDofs        :local DOFs for the velocity (9/27 per cell) 
%            - numHalfDofs     :the total number of HalfDofs (9/27*G.cells.num)
%            - HalfDofsToVdofs :an array where the number j in each row i
%                               corresponds to the Vdof of that HalfDof.

% $Date: 2009-10-01 17:06:53 +0200 (to, 01 okt 2009) $
% $Revision: 2929 $

   error(nargchk(1, 1, nargin, 'struct'));
   
   dim = size(G.nodes.coords,2);

   %{
  Nx = G.cartDims(1); Ny = G.cartDims(2);
  if dim == 2,
     [I, J   ] = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap);
  elseif dim == 3,
     Nz = G.cartDims(3);
     [I, J, K] = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap);
  else
     error(msgid('Dim:NotSupported'), ...
           '%d space dimensions not supported at this time.', dim);
  end
  if dim==3
    Nz = G.cartDims(3);
    [x{1:3}] = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap);
    G.cells.ijkMap = double([x{:}]);
    I = G.cells.ijkMap(:,1); J=G.cells.ijkMap(:,2); K=G.cells.ijkMap(:,3);
  else
    [x{1:2}] = ind2sub(reshape(G.cartDims, 1, []), G.cells.indexMap);
    G.cells.ijkMap = double([x{:}]);
    I=G.cells.ijkMap(:,1); J=G.cells.ijkMap(:,2);
  end
   %}

   % Numbering of nodes

   if dim == 2,
      % Numbering of nodes (4 per cell)
      % [lower left, lower right, upper left, upper right]
      nodes = cellNodes(G);
      nodeNums = reshape(nodes(:,3),4,[])';

      % Numbering of midpoints in each cell
      cellNums = (1:G.cells.num)';

      % Numbering of edge centers  - same as face numbering (4 per cell)
      % [left, right, lower, upper]
      edgeNums  = zeros(G.cells.num,4);
      edgeNums(:,1) = G.cellFaces((G.cellFaces(:,2)==1),1);
      edgeNums(:,2) = G.cellFaces((G.cellFaces(:,2)==2),1);
      edgeNums(:,3) = G.cellFaces((G.cellFaces(:,2)==3),1);
      edgeNums(:,4) = G.cellFaces((G.cellFaces(:,2)==4),1);

      Vdofs1 = nodeNums;                   % Nodal DOFs
      Vdofs2 = cellNums+max(max(Vdofs1));  % Cell DOFs
      Vdofs3 = edgeNums+max(Vdofs2);       % Edge center DOFs

      % Numbering of each velocity DOF per cell on each row
      Vdofs= [Vdofs1(:,1), Vdofs3(:,3), Vdofs1(:,2),...
              Vdofs3(:,1), Vdofs2(:),   Vdofs3(:,2),...
              Vdofs1(:,3), Vdofs3(:,4), Vdofs1(:,4)];

      % Pressure DOFs = the nodal DOFs
      Dofs.Pdofs    = Vdofs(:,[1,3,7,9]);
      Dofs.numPdofs = numel(unique(Dofs.Pdofs));

   elseif dim == 3,
      % Numbering of nodes
      %(8 per cell) [bottom lower left, bottom lower
      % right, bottom upper left, bottom upper right, top lower left, top
      % lower right, top upper left, top upper right]
      nodes = cellNodes(G);
      nodeNums = reshape(nodes(:,3),8,[])';

      % Numbering of midpoints in each cell
      cellNums = (1:G.cells.num)';

      % Numbering of edge centers - same as face numbering (6 per cell)
      % [left, right, back, front, bottom, top]
      faceNums      = zeros(G.cells.num,6);
      faceNums(:,1) = G.cellFaces((G.cellFaces(:,2)==1),1);
      faceNums(:,2) = G.cellFaces((G.cellFaces(:,2)==2),1);
      faceNums(:,3) = G.cellFaces((G.cellFaces(:,2)==3),1);
      faceNums(:,4) = G.cellFaces((G.cellFaces(:,2)==4),1);
      faceNums(:,5) = G.cellFaces((G.cellFaces(:,2)==5),1);
      faceNums(:,6) = G.cellFaces((G.cellFaces(:,2)==6),1);

      % Numbering of edge centers (12 per cell)
      %    [corners:l le, l r, u le, u r,...
      %             b lower, b upper, t lower, t upper ,...
      %             b left, b right, t left, t right]
      %
      edges = cellEdges(G);
      edgeNums       = zeros(G.cells.num,12);
      edgeNums(:, 1) = edges( 1:12:end,3);
      edgeNums(:, 2) = edges( 2:12:end,3);
      edgeNums(:, 3) = edges( 3:12:end,3);
      edgeNums(:, 4) = edges( 4:12:end,3);
      edgeNums(:, 5) = edges( 5:12:end,3);
      edgeNums(:, 6) = edges( 6:12:end,3);
      edgeNums(:, 7) = edges( 7:12:end,3);
      edgeNums(:, 8) = edges( 8:12:end,3);
      edgeNums(:, 9) = edges( 9:12:end,3);
      edgeNums(:,10) = edges(10:12:end,3);
      edgeNums(:,11) = edges(11:12:end,3);
      edgeNums(:,12) = edges(12:12:end,3);

      Vdofs1 = nodeNums;                   % Nodal DOFs
      Vdofs2 = cellNums+max(max(Vdofs1));  % Cell DOFs
      Vdofs3 = faceNums+max(max(Vdofs2));  % Face center DOFs
      Vdofs4 = edgeNums+max(max(Vdofs3));  % Edge center DOFs

      % Numbering of each velocity DOF per cell on each row
      Vdofs = [Vdofs1(:,1), Vdofs4(:,5), Vdofs1(:,2),...
               Vdofs4(:,9), Vdofs3(:,5), Vdofs4(:,10),...
               Vdofs1(:,3), Vdofs4(:,6), Vdofs1(:,4),...
               Vdofs4(:,1), Vdofs3(:,3), Vdofs4(:,2),...
               Vdofs3(:,1), Vdofs2(:),   Vdofs3(:,2),...
               Vdofs4(:,3), Vdofs3(:,4), Vdofs4(:,4),...
               Vdofs1(:,5), Vdofs4(:,7), Vdofs1(:,6),...
               Vdofs4(:,11),Vdofs3(:,6), Vdofs4(:,12),...
               Vdofs1(:,7), Vdofs4(:,8), Vdofs1(:,8)];

      % Pressure DOFs = the nodal DOFs
      Dofs.Pdofs    = Vdofs(:,[1,3,7,9,19,21,25,27]);
      Dofs.numPdofs = numel(unique(Dofs.Pdofs));
   end

  
  % Velocity DOFs
  Dofs.Vdofs    = Vdofs;
  Dofs.numVdofs = numel(unique(Dofs.Vdofs));
  
  if dim==2
    % HalfDofs (Velocity DOFs locally for each cell)
    HalfDofs         = reshape(1:G.cells.num*9,9,G.cells.num)';
    Dofs.HalfDofs    = HalfDofs;
    Dofs.numHalfDofs = G.cells.num*9;
    
    % PhalfDofs (Pressure DOFs locally for each cell)
    PhalfDofs         = reshape(1:G.cells.num*4,4,G.cells.num)';
    Dofs.PhalfDofs    = PhalfDofs;
    Dofs.numPhalfDofs = G.cells.num*4;
  elseif dim==3
    % HalfDofs (Velocity DOFs locally for each cell)
    HalfDofs         = reshape(1:G.cells.num*27,27,G.cells.num)';
    Dofs.HalfDofs    = HalfDofs;
    Dofs.numHalfDofs = G.cells.num*27;
    
    % PhalfDofs (Pressure DOFs locally for each cell)
    PhalfDofs         = reshape(1:G.cells.num*8,8,G.cells.num)';
    Dofs.PhalfDofs    = PhalfDofs;
    Dofs.numPhalfDofs = G.cells.num*8;
  end
  
  % HalfDofs2Vdofs
  HalfDofsToVdofs      = reshape(Dofs.Vdofs',[],1);
  Dofs.HalfDofsToVdofs = HalfDofsToVdofs;
  
  % PhalfDofs2Pdofs
  PhalfDofsToPdofs      = reshape(Dofs.Pdofs',[],1);
  Dofs.PhalfDofsToPdofs = PhalfDofsToPdofs;
  
  if dim==2
    % Face2Dofs
    cellFace2nodes      = reshape(Dofs.Vdofs(:,[1,4,7,3,6,9,1,2,3,7,8,9])',3,[])';
    Dofs.cellFace2nodes = cellFace2nodes;

    face2nodes                     = zeros(G.faces.num,3);
    face2nodes(G.cellFaces(:,1),:) = cellFace2nodes;
    Dofs.face2nodes               = face2nodes;
  elseif dim==3
    % Face2Dofs
    cellFace2nodes = reshape(Dofs.Vdofs(:,[ 1, 4, 7,10,13,16,19,22,25,...
                                            3, 6, 9,12,15,18,21,24,27,...
                                            1, 2, 3,10,11,12,19,20,21,...
                                            7, 8, 9,16,17,18,25,26,27,...
                                            1, 2, 3, 4, 5, 6, 7, 8, 9,...
                                           19,20,21,22,23,24,25,26,27])',9,[])';
    Dofs.cellFace2nodes = cellFace2nodes;
    
    face2nodes                     = zeros(G.faces.num,9);
    face2nodes(G.cellFaces(:,1),:) = cellFace2nodes;
    Dofs.face2nodes                = face2nodes;
  end
end
