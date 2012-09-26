function plotFaults(grdecl,G,varargin)
%  plot faults of the grid 
%
%  SYNOPSIS
%
%  plotFaults(grdecl,G,varargin)
%  
%  PARAMETERS:
%   grdecl - a grdecl structure
%   grid   - grid
%
%  RETURN:
%   none
% 
prm  = struct('radius',    5,  ...
              'height',    5,  ...
              'color',   'r',  ...
              'cylpts',   10,  ...
              'fontsize', 16,  ...
              'ambstr',  0.8);
prm = merge_options(prm, varargin{:});
faces=[];
color = colormap;
if(~isempty(grdecl))
    field_names = fieldnames(grdecl.FAULTS);
    for i = 1:length(field_names)
        field = char(field_names(i));
        % field is the fault tag
        faces = facesOfFault(field,grdecl,G);
        if(isfield(grdecl.FAULTS.(field),'multflt') && ~isempty(faces))
            multflt = grdecl.FAULTS.(field).multflt;
            multflt = min(multflt,1);
            name = field;
            plotFaces(G,faces,color(min(floor(multflt*64)+1,64),:));
            [val,ci] = min(G.faces.centroids(faces,3));
            pos = G.faces.centroids(ci,:);
            text(mean(G.faces.centroids(faces,1)),...
                mean(G.faces.centroids(faces,2)),...
                mean(G.faces.centroids(faces,3))-prm.height,...                 ...
                name, 'FontSize', prm.fontsize,  ...
                'Color', prm.color,'interp','none');
                   %'interp', 'latex',...
            %text(pos(1),pos(2),pos(3),name)
        end
    end
    colorbar
end