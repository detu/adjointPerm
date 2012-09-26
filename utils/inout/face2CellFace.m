function cellfaces = face2CellFace(faces,G)
%  calculate cell faces of given faces 
%
%  SYNOPSIS
%
%   cellfaces = face2CellFace(faces,G)
%  
%  PARAMETERS:
%   faces - list of faces
%   G     - grid struct
%
%  RETURN:
%   rock - stucture
%
%NB This is very slow should only be used if you really need this mapping
%  
disp('Warning: Very slow should not be used many faces')
cellfaces = zeros(length(faces),2);
for i=1:length(faces)
    ci=find(faces(i) == G.cellFaces(:,1));
    if(length(ci)==2)
        cellfaces(i,:) = ci;
    else
        cellfaces(i,:) = [ci, ci];
    end
end