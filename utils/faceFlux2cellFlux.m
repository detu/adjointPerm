function cellFlux = faceFlux2cellFlux(G, faceFlux)
%Transform face-based flux field to cell-based.
%
% SYNOPSIS:
%   cellFlux = faceFlux2cellFlux(G, faceFlux);
%
% PARAMETERS:
%   G        - Grid structure.
%
%   faceFlux - Vector of fluxes corresponding to face ordering.
%
% RETURNS:
%   cellFlux - Vector of fluxes cooresponding to cell face - ordering
%              corresponding to input.
%
% SEE ALSO:
%   cellFlux2faceFlux.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-14 13:22:35 +0200 (ma, 14 sep 2009) $
% $Revision: 2734 $

cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
sgn      = 2*(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo)-1;
cellFlux = sgn .* faceFlux(G.cellFaces(:,1));
