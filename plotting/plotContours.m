function plotContours(g, value, n)
% Plot contours of cell data.
%
% SYNOPSIS:
%   plotContours(g, value, n)
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
% RETURNS:
%
% NOTE:
%
% EXAMPLE
% 
% SEE ALSO:
%   plotCellData

%{
#COPYRIGHT#
%}

% $Id: $

levels = min(value): (max(value)-min(value))/n:max(value);
colors = round(1:size(colormap, 1)/n:size(colormap, 1));

for i=1:n,
   
   % fill space between contours
   I = [false; value >= levels(i) & value <= levels(i+1)];
   J = any(g.faces.neighbors==0, 2);
   f = find(any(I(g.faces.neighbors+1), 2) & J);
   plotFaces(g, f, colors(i), 'EdgeColor', 'none', 'facea', 0.4, 'outline', true);
   
   % Plot internal level sets
   I = [false; value <levels(i+1)];
   J = [false; value >levels(i+1)];
   f = find(any(I(g.faces.neighbors+1), 2) & any(J(g.faces.neighbors+1), 2));
   plotFaces(g, f, colors(i), 'EdgeColor', 'none', 'facea', 0.4, 'outline', true);
end

 
