function g = simpleGrdecl(dims, drop)
%Make a GRDECL structure for simple corner-point grid, possibly faulted.
%
% SYNOPSIS:
%   grdecl = simpleGrdecl(n)
%   grdecl = simpleGrdecl(n, drop)
%
% DESCRIPTION:
%   This function creates a nx-by-ny-by-nz cell corner-point grid.  The
%   resulting geometry discretises the physical domain
%
%      [0,1] x [0,1] x [0,0.5].
%
%   which is deformed by adding sin(x) and sin(x*y) perturbations.  In
%   addition, the pillars are skewed.  If a second input parameter is
%   provided, a single fault of given drop size is added in the model at
%   approximately x=0.5.
%
% PARAMETERS:
%   n    - Three-element vector, [nx, ny, nz], specifying the number of
%          cells in the 'x', 'y', and 'z' coordinate directions
%          respectively.
%
%   drop - Length of fault drop (i.e., the physical length (in metres) by
%          which the two blocks should be offset (along the fault)).
%          The drop size can be a numeric value or a handle to a function
%          of a spatial coordinate along the fault.
%          OPTIONAL.  Default value: drop = 0 (no fault).
%
% RETURNS:
%   grdecl - A GRDECL structure suitable for further processing by
%            function 'processGRDECL' or for permanent storage on disk by
%            function 'writeGRDECL'.
%
% EXAMPLE:
%   % Create a 20-by-20-by-5 fault grid file with a fault drop of 0.15 (m)
%   %grdecl = simpleGrdecl([20, 20, 5], @(x) 0.05*(sin(2*pi*x)-0.5) );
%   grdecl = simpleGrdecl([20, 20, 5], 0.15);
%
%   % Create the grid data structure
%   G = computeGeometry(processGRDECL(grdecl));
%
%   % Plot the geometry
%   hg = plotGrid(G, 'FaceAlpha', 0.625);
%   view(3), grid on, axis tight
%
% SEE ALSO:
%   processGRDECL, writeGRDECL.

%{
#COPYRIGHT#
%}

% $Id: simpleGrdecl.m 2596 2009-08-25 10:41:18Z bska $

%% Preprocess input data to make the underlying logical Cartesian grid

physDims   = [1, 1, 0.5];
g.cartDims = reshape(dims, 1, []);
[X, Y, Z]  = ndgrid(linspace(0, physDims(1), dims(1) + 1), ...
                    linspace(0, physDims(2), dims(2) + 1), ...
                    linspace(0, physDims(3), dims(3) + 1));

% Make skew pillars and perturb the layers.
X = X + 0.2*(0.5 - abs(X - 0.5)).*(Z - 0.5);

fun = @(x,y) -0.05*sin(pi*(x-0.5)) -0.075*sin(pi*(y+2*x));

for k = 1 : dims(3) + 1,
   xi       = X(:,:,k) ./ physDims(1);
   eta      = Y(:,:,k) ./ physDims(2);
   Z(:,:,k) = Z(:,:,k) + fun(xi, eta);
end
Z = sort(Z, 3);
Z = Z .* physDims(3);

%% Make pillars
n = prod(dims(1:2) + 1);
lines = zeros([n, 6]);
lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
g.COORD = reshape(lines.', [], 1);

%% Assign z-coordinates
% ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
z   = Z(ind(1), ind(2), ind(3));

%% Add fault
if nargin > 1,
   if isnumeric(drop),
      z(end/2+1:end,:,:) = z(end/2+1:end,:,:) + drop;
   elseif isa(drop, 'function_handle'),
      z(end/2+1:end,[1:2:end,end],:) = ...
         bsxfun(@plus, z(end/2+1:end,[1:2:end,end],:) ,...
                       drop(linspace(0,physDims(2), dims(2)+1)));
      z(end/2+1:end,2:2:end-1,:) =z(end/2+1:end,3:2:end,:);
   else
      error('');
   end
end
g.ZCORN = z(:);

%% Assign active cells

g.ACTNUM = reshape(ones(dims, 'int32'), [], 1);
