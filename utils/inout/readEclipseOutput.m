function output = readEclipseOutput(filename, varargin)
%read formated outputfiles from eclipse F* files
%
% SYNOPSIS:
% output = readEclipseOutput(filename)
% 
% PARAMETERS:
%   filename - string with the filename of the file
%
% RETURNS:
%   output - struct with all the data  in the file
%
% SEE ALSO:
%   readEclipseOutputField.m,readEclipseTimeData.m,readEclipseResult.m
%
clear output
fid = fopen(filename,'r');
if(fid == -1)
    error(['File: ',filename,' do not exist']);
end
while ~feof(fid),
    [name,field]=readEclipseOutputField(fid, varargin{:});
    ind = find(name == '-');
    if(~isempty(ind))
        name(ind) = '_';
    end
    if(~isempty(field))        
        output.(name)=field;
    end    
end
fclose(fid);
end
