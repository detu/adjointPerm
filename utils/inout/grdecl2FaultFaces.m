function faces = grdecl2FaultFaces(grdecl,G)
faces=[];
if(~isempty(grdecl))
    field_names = fieldnames(grdecl.FAULTS);
    for i = 1:length(field_names)
        field = char(field_names(i));
        % field is the fault tag
        faces = [faces;facesOfFault(field,grdecl,G)];
    end
end