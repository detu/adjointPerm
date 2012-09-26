function [well_data,cell_data] = readEclipseResults(file_prefix,varargin)
% read formatted output from eclipse directory and output structs for
% welldata and celldata as time series
%
% SYNOPSIS:
%  [well_data,cell_data] = readEclipseResults(file_prefix)
%  [well_data,cell_data] = readEclipseResults(file_prefix,varargin)
% PARAMETERS:
%   file_prefix - name of file prefix time files of  file_prefix.A0* for
%                 production data and file_prefix.F* for celldata should
%                 exist
% OPTIONAL                  
%   files_out - give file name instead of data, default is true (to save
%               space
%
% RETURNS:
%   well_data - struct with the given welldata defined
%   cell_data - structure of celldata with the well_data.celldata{i} as
%               as the celldata at time cell_data.times(i)
% SEE ALSO:
%   readEclipseOutputField.m,readEclipseOutput.m,readEclipseTimeData.m
%
%file_prefix = 'new_data/BC0407_IO'
opt = struct('verbose' , false, ...
             'files_out',true);             
opt = merge_options(opt, varargin{:});
dir_prefix = regexp(file_prefix, '[''"]?([\w^\d]+)[''"]?','tokens', 'once');
file = [file_prefix,'.FSMSPEC'];
well_def =readEclipseOutput(file);

simfiles = dir([file_prefix,'.A*']);
cellfiles = dir([file_prefix,'.F0*']);
if(isempty(simfiles))
    error('No files')
end
simnames = {simfiles.name};

%fix sorting of files not konsitent with time ??
%pattern = '^\w+_(\d+)_(\d+)_(\d+)';
%N = regexp(simnames, pattern, 'tokens');
%[sort_ind,sort_ind] = sortrows(cell2mat(cellfun(@str2double, cat(1, N{:}), 'unif', false)));
%data1 =readEclipseOutput(filename)
%num_files= numel(simnames);
sort_ind = 1:numel(simnames);

%cell_data = struct('times',-ones(num_files,1),{},'saturation',{})
well_data=[];
cell_data =[];
for i=sort_ind
    filename = [char(dir_prefix),'/',char(simnames(i))];
    welldata = readEclipseTimeData(filename,well_def);
    filename = [char(dir_prefix),'/',cellfiles(i).name];
    if(isempty(well_data))
        well_data = welldata;
        cell_data.times = welldata.TIME(end);
    else
        field_names = fieldnames(welldata);
        for j = 1:numel(field_names);
            field = char(field_names(j));
            if(~strcmp(field,'WNAME'))
                well_data.(field) = [well_data.(field),welldata.(field)];
            end
        end
        cell_data.times = [cell_data.times,welldata.TIME(end)];
    end
    if(~opt.files_out)
        celldata = readEclipseOutput(filename);
        cell_data = addCellData(cell_data,celldata);
    end
    if(isfield(cell_data,'files'))
        cell_data.files{end+1} = filename;
    else
        cell_data.files ={filename};
    end
end
end
function cell_data = addCellData(cell_data,celldata)
   %fieldnames = fields(celldata); 
   fieldnames = {'PRESSURE','SWAT'};
   for i=1:numel(fieldnames)
       name = char(fieldnames(i));
       if(isfield(cell_data,fieldnames(i)))
        cell_data.(name){end+1} = celldata.(name).values;
       else
          cell_data.(name) = {celldata.(name).values};
       end
   end
end
