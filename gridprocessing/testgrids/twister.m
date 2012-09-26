function pt = twister(pt) 
%Permutes x- and y-coordinates of nodes in a grid.
%
% SYNOPSIS:
%  twister(G.nodes.coords)
%
% PARAMETERS:
%  pt - Coordinates to permute
% 
% RETURNS:
%  pt - Modified coordinates

%{
#COPYRIGHT#
%}

%$Id:$

xi = bsxfun(@rdivide, bsxfun(@minus, pt, min(pt)), max(pt)-min(pt));
f  = @(x,y)0.03*sin(pi*x).*sin(3*(-pi/2+pi*y));
pt(:,1:2) = bsxfun(@times, [xi(:,1)+f(xi(:,1), xi(:,2)),...
                   xi(:,2)-f(xi(:,2), xi(:,1))], ...
                   max(pt(:,1:2))-min(pt(:,1:2)));
end
