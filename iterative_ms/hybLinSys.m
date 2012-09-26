function [B, P, f, g, h, mobt] = ...
    hybLinSys(resSol, G, S, rock, fluid, p0, dt, bc, src)


%--------------------------------------------------------------------------
% Get updated pressure and saturation dependent quantities -------------
[c, rho, mu, u] = fluid.pvt(resSol.cellPressure, resSol.z);
s      = bsxfun(@rdivide, u, sum(u,2));
ct     = sum( s  .*  c , 2);
mob    = fluid.relperm(s) ./ mu;
mobt   = sum(mob, 2);
Lti    = 1./mobt;
rhoTot = sum(rho.*mob ,2).*Lti;

%--------------------------------------------------------------------------
% System matrices ---------------------------------------------------------
nc     = G.cells.num;
cf     = G.cellFaces(:, 1);
ncf    = numel(cf);
cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';

B       = sparse(1:ncf, 1:ncf, Lti(cellNo)) * S.B;
accum   = ct .* poreVolume(G, rock) ./ dt;
P       = sparse(1:nc,1:nc,-accum);

%--------------------------------------------------------------------------
% System RHS --------------------------------------------------------------
grav = gravity(); grav = grav(:);

[f, g, h, fg] = computePressureRHS(G, rhoTot, bc, src);
f = f + fg;

% Accumulation
g = g + accum .* p0;

% a1 - part
a1 = sum( c .* mob, 2) .* Lti;
%vv = sparse((1:ncf)',cellNo, resSol.cellFlux, ncf, nc);
%g  = g + a1 .* ( vv' * ( B*resSol.cellFlux ) );
g = g + a1.*accumarray(cellNo, resSol.cellFlux.*(B*resSol.cellFlux), [nc,1]);

if any(grav)
   % a2 - part
   hg = ( G.faces.centroids(cf, :) - G.cells.centroids(cellNo, :) ) * grav;
   a2 = sum( c  .* mob .* bsxfun(@minus, 2*rhoTot, rho), 2) .* Lti;  % sjekk fortegn???
   % g  = g - a2 .* (vv'*hg);
   g  = g - a2 .* accumarray(cellNo,resSol.cellFlux.*hg, [nc,1]);
   % a3 - part
   a3 = rhoTot.* sum( c .* mob .* bsxfun(@minus, rhoTot, rho));     % sjekk fortegn???
   [K, row, col] = permTensor(rock, size(G.nodes.coords,2));
   gKg  = bsxfun(@times, K, grav(col)) .* grav(row)';
   g = g + a3 .* poreVolume(G, rock) .* gKg;
end
