function [htop, htext, hs] = plotWell(G, W, varargin)
%Plot well trajectories into current axes.
%
% SYNOPSIS:
%   [htop, htext, hs] = plotWell(G, W)
%   [htop, htext, hs] = plotWell(G, W, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure
%   W       - Well data structure
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%             Supported options are:
%               'radius' -- Real, > 0, denoting the width of the well path
%                           when plotting into the axes.
%                           Default value: radius = 5
%               'height' -- Height above top reservoir contact at which the
%                           well should stop and symbold drawn.
%                           Default value: height = 5
%               'color'  -- Colour with which the well path should be
%                           drawn.  Possible values described in PLOT.
%                           Default value: color = 'r'
%               'cylpts' -- Number of segments to use about the well bore.
%                           A higher value produces more smoothly looking
%                           well trajectories at the expense of more costly
%                           plotting.
%                           Default value: cylpts = 10
%             'fontsize' -- The size of the font used for the label texts.
%                           Default value: fontsize = 16
%               'ambstr' -- The ambient strenght of the well cylinder.
%                           Default value: ambstr = 0.8
%
% RETURNS:
%   htop  - Graphics handles to all well heads and heels.
%   htext - Graphics handles to all rendered well names.
%   hs    - Graphics handles to all well paths.
%
% EXAMPLE:
%   See simpleWellExample.m.
%
% SEE ALSO:
%   addWell, delete, patch.

%{
#COPYRIGHT#
%}

% $Id: plotWell.m 2637 2009-09-03 09:27:13Z hnil $

prm  = struct('radius',    5,  ...
              'height',    5,  ...
              'color',   'r',  ...
              'cylpts',   10,  ...
              'fontsize', 16,  ...
              'ambstr',  0.8);
prm = merge_options(prm, varargin{:});

dohold = ishold;
hold on;

nW    = numel(W);
% Number of digits required to textually represent all wells on the form
% W%d (e.g. 'W03') with the same number of digits for all wells.
%
nWd = fix(log10(nW)) + 1;

htop  = zeros([1, 2*nW]);
htext = zeros([1,   nW]);
hs    = zeros([1,   nW]);

% Common face/vertex connection table for plotting well paths by means of
% PATCH.
faces = @(sz,i,j) [sub2ind(sz, i  , j  ), ... % lower left  vertex
                   sub2ind(sz, i+1, j  ), ... % lower right vertex
                   sub2ind(sz, i+1, j+1), ... % upper right vertex
                   sub2ind(sz, i  , j+1)];    % upper left  vertex

for w = 1 : nW,
   % Plot cylinder
   c = G.cells.centroids(W(w).cells, :);

   dir = W(w).dir;
   if numel(W(w).cells) == 1,
      c = c([1 1],:);
      dir = 's';
   end

   % Generate a cylinder rotatet according to major direction
   switch dir
      case 'x',
         [zw, yw, xw] = cylinder(prm.radius(ones([size(c,1),1])), prm.cylpts);
         xw = 0*xw;
      case 'y',
         [xw, zw, yw] = cylinder(prm.radius(ones([size(c,1),1])), prm.cylpts);
         yw = 0*yw;
      case 'z',
         [xw, yw, zw] = cylinder(prm.radius(ones([size(c,1),1])), prm.cylpts);
         zw = 0*zw;
      case 's',
         [xw, yw, zw] = cylinder(prm.radius(ones([size(c,1),1])), prm.cylpts);
         zw = 4*prm.radius*(zw-0.5);
       otherwise
          [xw, yw, zw] = cylinder(prm.radius(ones([size(c,1),1])), prm.cylpts);
          zw = 4*prm.radius*(zw-0.5);
   end

   xw      = xw + repmat(c(:,1), [1, prm.cylpts + 1]);
   yw      = yw + repmat(c(:,2), [1, prm.cylpts + 1]);
   zw      = zw + repmat(c(:,3), [1, prm.cylpts + 1]);
   %zw(1,:) = zw(1,:) - prm.height;

   % Draw top
   nodes = [xw(1,:); yw(1,:); zw(1,:)] .';
   htop(2*(w-1) + 1) = patch('Faces',  (1 : numel(nodes)/3), ...
                             'Vertices',  nodes,             ...
                             'EdgeColor', 'k',               ...
                             'FaceColor', prm.color,         ...
                             'AmbientStrength', prm.ambstr);

   % Draw bottom
   nodes = [xw(end,:); yw(end,:); zw(end,:)] .';
   htop(2*w) = patch('Faces',  (1 : numel(nodes)/3), ...
                     'Vertices', nodes,              ...
                     'EdgeColor', 'k',               ...
                     'FaceColor', prm.color,         ...
                     'AmbientStrength', prm.ambstr);

   if isfield(W(w), 'name') && ischar(W(w).name),
      name = W(w).name;
   else
      name = sprintf('W$_{%0*d}$', nWd, w);
   end
   htext(w) = text(mean(xw(1,:)),                   ...
                   mean(yw(1,:)),                   ...
                   mean(zw(1,:))-prm.height,        ...
                   name, 'FontSize', prm.fontsize,  ...
                   'interp', 'latex', 'Color', prm.color);

   sz    = size(xw);
   [I,J] = ndgrid(1 : sz(1)-1, 1 : sz(2)-1);
   hs(w) = patch('Faces'          , faces(sz, I(:), J(:)), ...
                 'Vertices'       , [xw(:), yw(:), zw(:)], ...
                 'AmbientStrength', prm.ambstr,            ...
                 'FaceColor'      , prm.color,             ...
                 'EdgeColor'      , 'none');

   % Plot a line to somewhere above the model
   plot3(c([1,1],1), c([1,1],2), c([1,1],3) - [0; prm.height], ...
         'color', prm.color, 'LineWidth', 2);
end

if dohold,
   hold on;
else
   hold off;
end
