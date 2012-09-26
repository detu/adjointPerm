function LS = initMixedBO(G, bc, src, W, flag) %#ok

LS.bc = bc;
LS.src = src;
LS.wells = W;

nc  = G.cells.num;
nf  = G.faces.num;
cf  = G.cellFaces(:, 1);

% Global to local faces
ncf = numel(cf);
cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';
sgn    = 2*( G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo ) -1;

C   = sparse((1:ncf)', cellNo, 1, ncf, nc);
D   = sparse((1:ncf)', double(cf), 1, ncf, nf);

if isempty(W)
   nperf = 0;
   Dw = deal([]);
   Cw = sparse(0, G.cells.num);
else
   nperf = cellfun(@numel, { W.cells });
   wc = vertcat(W.cells);
   np = numel(wc);
   nw = numel(W);
   Cw = sparse(1:np, wc , 1, np, G.cells.num);
   Dw = sparse(1:np, rldecode(1:nw, nperf, 2), 1, np, nw);
end

% Mapping: faceFlux -> cellFlux
LS.Do = blkdiag( sparse(1:ncf,1:ncf, sgn)*D, ...
              speye(sum(nperf), sum(nperf)) );
           
% Reduce to Neumann faces (D and h)
[fluxFacesR, fluxFacesW, pFacesR, pFacesW] = getBCType(G, W, bc);
LS.DN  = blkdiag( D(:, fluxFacesR), Dw(:, fluxFacesW) );
LS.fluxFacesR = fluxFacesR;
LS.fluxFacesW = fluxFacesW;
LS.pFacesR    = pFacesR;
LS.pFacesW    = pFacesW;

LS.C = [C; Cw];
if nargin>4, LS.D = blkdiag(D, Dw); end

LS.solved = true;  % needs to assemble dynamic parts before next solve