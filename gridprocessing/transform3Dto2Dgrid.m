function G2 = transform3Dto2Dgrid(G)
% Transforms a 3D grid into a 2D grid. 
% 3D grid must be constant in either x-, y- or z-direction.
  
  if (numel(unique(G.cells.centroids(:,1)))==1) 
    dir=1;
  elseif (numel(unique(G.cells.centroids(:,2)))==1) 
    dir=2;
  elseif (numel(unique(G.cells.centroids(:,3)))==1) 
    dir=3;
  end
  
  % G2.cells
  G2.cells.num         = G.cells.num;
  G2.cells.numFaces    = G.cells.numFaces;
  G2.cells.numFaces(:) = 4;
  if(isfield(G.cells,'indexMap'))	          
  	G2.cells.indexMap = G.cells.indexMap;
  end	
  
  if dir==1
    dx = max(G.nodes.coords(:,1))-min(G.nodes.coords(:,1));
    G2.cells.volumes   = G.cells.volumes./dx;
    G2.cells.centroids = G.cells.centroids(:,2:3);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,2))),...
                           numel(unique(G.cells.centroids(:,3)))]);
    
    % G2.nodes
    nodes           = find(G.nodes.coords(:,1)== min(G.nodes.coords(:,1)));
    G2.nodes.coords = G.nodes.coords(nodes,[2,3]); 
    G2.nodes.num    = size(G2.nodes.coords,1);

    % G2.faces
    faces = true(G.faces.num,1);
    faces(find(G.faces.normals(:,1))) = false;
    G2.faces.num         = numel(find(faces));
    G2.faces.numNodes    = G.faces.numNodes(find(faces));
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(find(faces),:);
    G2.faces.areas       = G.faces.areas(find(faces))./dx;
    G2.faces.normals     = G.faces.normals(find(faces),2:3);
    G2.faces.centroids   = G.faces.centroids(find(faces),2:3);
    
    % G2.faceNodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faceNodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faceNodes = facenodes(ismember(facenodes,nodes));
    G2.faceNodes = compressFaces(G2.faceNodes); 
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                             double(G2.faces.numNodes)), G2.faceNodes];
    y_faces      = find(G2.faces.normals(:,2)); 
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:);
    y_face       = 0;
    for i = 1:size(y_faces) 
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;  
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
     G2.faceNodes(ind_fn(y_face)) = G2.faceNodes(ind_fn(y_face(ind)));
    end
    y_faces      = find(G2.faces.normals(:,1));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,1));
    y_facenodes  = facenodes(ind_fn,:); 
    y_face       = 0;
    for i = 1:size(y_faces) 
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;  
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,2), 'ascend');
      G2.faceNodes(ind_fn(y_face)) = G2.faceNodes(ind_fn(y_face(ind)));
    end
    
    % G2.cellFaces
    faces        = true(size(G.cellFaces,1),1);
    faces(find(G.cellFaces(:,2)==1|G.cellFaces(:,2)==2)) = false;
    G2.cellFaces = G.cellFaces(find(faces),:);  
    G2.cellFaces = compressFaces(G2.cellFaces);  
    x_faces1 = find(G2.cellFaces(:,2)==3);
    x_faces2 = find(G2.cellFaces(:,2)==4);
    x_faces3 = find(G2.cellFaces(:,2)==5);
    x_faces4 = find(G2.cellFaces(:,2)==6);
    G2.cellFaces(x_faces1,2)=1;
    G2.cellFaces(x_faces2,2)=2;
    G2.cellFaces(x_faces3,2)=3;
    G2.cellFaces(x_faces4,2)=4;

  elseif dir==2
    dy = max(G.nodes.coords(:,2))-min(G.nodes.coords(:,2));
    G2.cells.volumes   = G.cells.volumes./dy;
    G2.cells.centroids = G.cells.centroids(:,[1,3]);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,1))),...
                           numel(unique(G.cells.centroids(:,3)))]);

    % G2.nodes
    nodes           = find(G.nodes.coords(:,2)== min(G.nodes.coords(:,2)));
    G2.nodes.coords = G.nodes.coords(nodes,[1,3]); 
    G2.nodes.num    = size(G2.nodes.coords,1);
    
    % G2.faces
    faces = true(G.faces.num,1);
    faces(find(G.faces.normals(:,2))) = false;
    G2.faces.num         = numel(find(faces));
    G2.faces.numNodes    = G.faces.numNodes(find(faces));
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(find(faces),:);
    G2.faces.areas       = G.faces.areas(find(faces))./dy;
    G2.faces.normals     = G.faces.normals(find(faces),[1,3]);
    G2.faces.centroids   = G.faces.centroids(find(faces),[1,3]);

    % G2.faceNodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faceNodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faceNodes = facenodes(ismember(facenodes,nodes));
    G2.faceNodes = compressFaces(G2.faceNodes); 
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                             double(G2.faces.numNodes)), G2.faceNodes];
    y_faces      = find(G2.faces.normals(:,2)); 
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:); 
    y_face       = 0;
    for i = 1:size(y_faces) 
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;  
      node       = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
      G2.faceNodes(ind_fn(y_face)) = G2.faceNodes(ind_fn(y_face(ind)));
    end

    % G2.cellFaces
    faces        = true(size(G.cellFaces,1),1);
    faces(find(G.cellFaces(:,2)==3|G.cellFaces(:,2)==4)) = false;
    G2.cellFaces = G.cellFaces(find(faces),:);  
    G2.cellFaces = compressFaces(G2.cellFaces);  
    y_faces1     = find(G2.cellFaces(:,2)==5);
    y_faces2     = find(G2.cellFaces(:,2)==6);
    G2.cellFaces(y_faces1,2) = 3;
    G2.cellFaces(y_faces2,2) = 4;
	
elseif dir==3
    dz = max(G.nodes.coords(:,3))-min(G.nodes.coords(:,3));

    G2.cells.volumes   = G.cells.volumes./dz;   
    G2.cells.centroids = G.cells.centroids(:,1:2);
    G2.cartDims        = ([numel(unique(G.cells.centroids(:,1))),...
                           numel(unique(G.cells.centroids(:,2)))]);

    % G2.nodes
    nodes           = find(G.nodes.coords(:,3)== min(G.nodes.coords(:,3)));
    G2.nodes.coords = G.nodes.coords(nodes,1:2); 
    G2.nodes.num    = size(G2.nodes.coords,1);

    % G2.faces
    faces = true(G.faces.num,1);
    faces(find(G.faces.normals(:,3))) = false;
    G2.faces.num         = numel(find(faces));
    G2.faces.numNodes    = G.faces.numNodes(find(faces));
    G2.faces.numNodes(:) = 2;
    G2.faces.neighbors   = G.faces.neighbors(find(faces),:);
    G2.faces.areas       = G.faces.areas(find(faces))./dz;
    G2.faces.normals     = G.faces.normals(find(faces),1:2);
    G2.faces.centroids   = G.faces.centroids(find(faces),1:2);

    % G2.faceNodes:
    facenodes    = [rldecode((1 : G.faces.num) .', ...
                             double(G.faces.numNodes)), G.faceNodes];
    actfaces     = ismember(facenodes(:,1),find(faces));
    facenodes    = facenodes(actfaces,2);
    G2.faceNodes = facenodes(ismember(facenodes,nodes));
    facenodes    = [rldecode((1 : G2.faces.num) .', ...
                                     double(G2.faces.numNodes)), G2.faceNodes];
    y_faces      = find(G2.faces.normals(:,2));
    norm_fn      = rldecode(G2.faces.normals,double(G2.faces.numNodes));
    ind_fn       = find(norm_fn(:,2));
    y_facenodes  = facenodes(ind_fn,:); 
    y_face       = 0;
    for i=1:size(y_faces)
      face       = y_faces(i);
      y_face     = max(y_face)+1:max(y_face)+1+...
                   double(G2.faces.numNodes(face))-1;  
      node = y_facenodes(y_face,2);
      [srt, ind] = sort(G2.nodes.coords(node,1), 'descend');
      G2.faceNodes(ind_fn(y_face)) = G2.faceNodes(ind_fn(y_face(ind)));
    end

    % G2.cellFaces
    faces        = true(size(G.cellFaces,1),1);
    faces(find(G.cellFaces(:,2)==5|G.cellFaces(:,2)==6)) = false;
    G2.cellFaces = G.cellFaces(find(faces),:);  
    G2.cellFaces = compressFaces(G2.cellFaces);  
  end
end
  
function cf = compressFaces(cf);
  f=cf(:,1);
  active = find(accumarray(f, 1) > 0);
  compr  = zeros([max(f), 1]);
  compr(active) = 1 : numel(active);
  f = compr(f);
  cf(:,1)=f;
end
