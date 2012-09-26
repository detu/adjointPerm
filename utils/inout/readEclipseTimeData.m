function welldata = readEclipseTimeData(filename,well_def)
%read formatted output of well data from eclise outputfile (A* files)
%
% SYNOPSIS:
%  welldata = readEclipseTimeData(filename,well_def)
% 
% PARAMETERS:
%   filename - name of file  (*A0*)
%   well_def - struct read from the a definition file 'name.FSMSPEC'
%     well_def =readEclipseOutput(file);
%
%
% RETURNS:
%   welldata - struct with the given welldata defined
%
% SEE ALSO:
%   readEclipseOutput.m,readEclipseOutputField.m,readEclipseResult.m
%
output=[];
clear output
fid = fopen(filename,'r');
if(fid == -1)
    error(['File: ',filename,' do not exist']);
end
count =1; 
output.MINISTEP =[];
while ~feof(fid),
    [name,field]=readEclipseOutputField(fid);    
    if(strcmp(name,'MINISTEP'))        
        number = field.values;
        [name,field]=readEclipseOutputField(fid) ;
        if(number == 0)
            output.(name) = {};
        end
        output.MINISTEP =[output.MINISTEP,number];
        output.(name){count} = field;
        count =count+1;
    else
        output.(name) = field;
    end
end
fclose(fid);
%%
numwells=well_def.DIMENS.values(4);
wellpos=[];
for i=1:numwells
    wellpos=[wellpos,find(well_def.NUMS.values == i)];
end


%
% read out not well spesific data
for j = 1:19
   name = well_def.KEYWORDS(j).values; 
    name = char(regexp(name,'["]?([-./\w]+)[''"]?','tokens','once'));
    welldata.(name) =[];
    for k=1:size(output.PARAMS,2)
      welldata.(name) = [welldata.(name),output.PARAMS{k}.values(j)];
    end
end

%% start reading well data
%find position of well data ? 
for i=1:numwells
    name = well_def.WGNAMES(wellpos(1,i)).values;
    name = char(regexp(name,'["]?([-./\w]+)[''"]?','tokens','once'));
    welldata.WNAME{i} = name;
end

% read well data
for j=1:size(wellpos,1)
    %well_def.KEYWORDS(reshape(wellpos(j,:),1,[])).values %(all equal ?    
    name = well_def.KEYWORDS(wellpos(j,1)).values;
    name = char(regexp(name,'["]?([-./\w]+)[''"]?','tokens','once'));
    welldata.(name) =[];
    for k=1:size(output.PARAMS,2)
         welldata.(name) = [welldata.(name),output.PARAMS{k}.values(wellpos(j,:))];
    end
    %welldata.(name) = welldata.(name)';
end
welldata.MINISTEP = output.MINISTEP;