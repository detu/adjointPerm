function [fluxFacesR, fluxFacesW, pFacesR, pFacesW] = getBCType(G, W, bc)
%
% Utility function for internal use only
%
fluxFacesR = any(G.faces.neighbors == 0, 2);
pFacesR    = false(size(fluxFacesR));
if ~isempty(bc)
    pf    = strcmpi('pressure', bc.type');
    pIndx = double(bc.face(pf));
    fluxFacesR ( pIndx ) = false;
    pFacesR( pIndx ) = true;
end
if nargout > 1
    if ~isempty(W)
        fluxFacesW = strcmpi('rate', {W.type}');
        pFacesW = strcmpi('pressure', {W.type}');
    else
        fluxFacesW = logical([]);
        pFacesW = logical([]);
    end
end