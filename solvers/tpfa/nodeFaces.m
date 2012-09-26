function G = nodeFaces(G)

   G.nodes.numFaces = accumarray(G.faceNodes, 1, [G.nodes.num, 1]);
   assert(all( G.nodes.numFaces > 0));
   
   faceNo = rldecode(1:G.faces.num, double(G.faces.numNodes), 2) .';
   
   tab = sortrows([G.faceNodes, faceNo]);
   G.nodeFaces = uint32(tab(:,2));
   
   
   
