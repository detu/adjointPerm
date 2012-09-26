function faceFlux = cellFlux2faceFlux(G, cellFlux)
%Transform cell-based flux field to face-based.
%
% SYNOPSIS:
%   faceFlux = cellFlux2faceFlux(G, cellFlux);
%
% PARAMETERS:
%   G        - grid structure
%
%   cellFlux - Vector of fluxes corresponding to cell-face ordering
%
% RETURNS:
%   faceFlux - Vector of fluxes corresponding to face-ordering
%              corresponding to input.
%
% SEE ALSO:
%   faceFlux2cellFlux.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-14 13:22:35 +0200 (ma, 14 sep 2009) $
% $Revision: 2734 $

cf       = G.cellFaces(:,1);
cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
sgn      = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;
faceFlux = accumarray(cf, sgn .* cellFlux, [G.faces.num, 1]) ./ ...
           accumarray(cf,       1        , [G.faces.num, 1]);

faceFlux(~isfinite(faceFlux)) = 0;
