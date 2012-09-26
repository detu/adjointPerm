function [name,output] = readEclipseOutputField(fid, varargin)
%read the next field of a formated eclipe output file (F* files)
%
% SYNOPSIS:
% [name,output] = readEclipseOutputField(fid)
% 
% PARAMETERS:
%   fid - open file id
%
%
% RETURNS:
%   output - struct with all the data  field
%   name   - of the field
%
% SEE ALSO:
%   readEclipseOutput.m,readEclipseTimeData.m,readEclipseResult.m
%

opt = struct('verbose', false);
opt = merge_options(opt, varargin{:});

lin = fgetl(fid);
a = regexp(lin, '[''"]?([\.\w\-/]+)[''"]?','tokens');
    name = a{1}{:};
    number = str2num(a{2}{:});
    ttype = a{3}{:};
if strcmp(ttype,'INTE'),
    values = textscan(fid,'%f',number);
elseif strcmp(ttype,'DOUB') || strcmp(ttype,'REAL'),
    values = textscan(fid,'%f',number);
elseif strcmp(ttype,'LOGI'),
    values = textscan(fid,'%s',number);
elseif strcmp(ttype,'LOGI') || strcmp(ttype,'CHAR'),
    %disp('char')
    values =textscan(fid,'%s',number,'Delimiter','''','MultipleDelimsAsOne',true);
else
    if number == 0,
        dispif(opt.verbose, 'Keyword ''%s'' with zero elements', name);
        output = struct('type',ttype,'values',{});
        return
    else
        error('Unknow variable');
    end
end
lin = fgetl(fid);%#ok
output = struct('type',ttype,'values',values{:});

