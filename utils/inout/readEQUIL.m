function grdecl = readEQUIL(fid,grdecl)
% Read EQUIL keyword
%
%  SYNOPSIS
%
%   grdecl = readEQUIL(fid,grdecl)
%  
%  PARAMETERS:
%   grdecl - a grdecl structure
%
%  RETURN:
%   grdecl - a grdecl structure
%
lin = fgetl(fid);
count = 1;
if(~isfield(grdecl,'EQUIL'))
    grdecl.EQUIL=[];
end
while isempty(regexp(lin, '^/', 'once'))
    % Skip blank lines and comments.
    if(~(isempty(lin) || all(isspace(lin)) || ~isempty(regexp(lin, '^--', 'match'))))
        split = regexp(lin, '(\d\.*\-*\**)+', 'match');
        if(~length(split)==9)
            error(['Wrong line in readOpertor ', kw ,' for ',name])
        end
        grdecl.EQUIL = [grdecl.EQUIL;str2double(split)];
    end
    lin = fgetl(fid);
    count =count +1; 
end
