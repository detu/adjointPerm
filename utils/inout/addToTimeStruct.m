function mt_struct = addToTimeStruct(mt_struct,t_struct)
% extend all stuct names with a stuct with the same fields defined
%
% SYNOPSIS:
% output = readEclipseOutput(filename)
% 
% PARAMETERS:
%   mt_struct - multi time stuct
%   t_struct - one time stuct
%
% RETURNS:
%   mt_struct - same stuct with extra added by [,]
%
field_names = fieldnames(t_struct);
for i=1:numel(field_names)
       name = char(field_names(i));
       if(~isfield(mt_struct,field_names(i)))
            mt_struct.(name) = t_struct.(name);
       else
          mt_struct.(name) = [mt_struct.(name),t_struct.(name)];
       end
   end
end