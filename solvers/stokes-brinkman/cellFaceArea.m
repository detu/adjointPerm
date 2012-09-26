function actArea = cellFaceArea(CS, G, CG)
% cellFaceArea -- Calculates the area of the active cell faces CS.activeCellFaces.
% 
% SYNOPSIS:
%   actArea = cellFaceArea(CS, G, CG)
%
% PARAMETERS:
%
%   CS    - Coarse system structure.
%
%   G, CG - Fine/coarse grid structure as generated.
%
%
% RETURNS:
%   actArea - Array of (CS.activeCellFaces,1) containing the area of the
%   active cell faces . 
 
    
  [nsub, sub] = subFaces(G, CG);
  sub_ix      = cumsum([0; nsub]);
  faces=[1:CG.faces.num];
  FaceArea=arrayfun(@(actface) sum(G.faces.areas(sub(sub_ix(actface) + 1:...
      sub_ix(actface + 1)))),faces);
  
  facenum=false(CG.faces.num,1);   
  facenum(CS.activeFaces)=true;    
 
  actArea=zeros(numel(CS.activeCellFaces),1);
  actArea=FaceArea(CG.cellFaces(facenum(CG.cellFaces(:,2)),2))';
end