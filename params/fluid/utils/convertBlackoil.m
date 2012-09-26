function [c, rho, u] = convertBlackoil(z, B, dBdp, R, dRdp, rho_S, col)
%Compute compressibility, density and phase volume from blackoil data
%
% SYNOPSIS:
%   [d, rho, u] = convertBlackoil(z, B, dBdp, R, dRdp, rho_S, col)
%
% PARAMETERS:
%   z       - Surface volumes
%   B       - Formation volume factors
%   dBdp    - B differentiated wrt pressure
%   R       - Miscibility ratios (gas-liquid, oil-vapour and gas-aqua)
%   dRdp    - R differentiated wrt pressure
%   rho_S   - Surface density
%   col     - Sequence of phases. {A,L,V} where A is column nnumber for the
%             aqueous phase, L is the column number for the liquid phase
%             and V is the column number for the vapour phase.
%
% RETURNS:
%   c       - Phase compressibilities
%   rho     - Phase densities
%   u       - Phase volumes
%
% EXAMPLE:
%
% SEE ALSO:
%   initBlackoilFluid, initBlackOilRelPerm.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-10 10:19:46 +0200 (on, 10 jun 2009) $
% $Revision: 2354 $

   c    = phaseCompressibility(B, dBdp, R, dRdp, col);
   rho  = phaseDensity(B, R, rho_S, col);
   u    = phaseVolume(z, B, R, col);
end


% R(:,L) = Rs = R_Liquid  = gas dissolved in oil
% R(:,V) = Rv = R_Vapor   = oil vaporized in gas
% R(:,A) = Ra = R_Aqua    = gas dissolved in water
function u   =  phaseVolume(z, B, R, col)
   %u = B*inv(R)*z;
   A = col{1};  L = col{2};  V = col{3};

   if A, u(:,A)   = B(:,A).* z(:,A);end
   if L, u(:,L)   = B(:,L).* z(:,L);end
   if V, u(:,V)   = B(:,V).* z(:,V);end

   % If we allow mixing, both solute and solvent must be present.
   if L & V, %#ok
      u(:,L) = u(:,L) - B(:,L).*R(:,V).*z(:,V);
      u(:,V) = u(:,V) - B(:,V).*R(:,L).*z(:,L);
   end

   % If we allow mixing, both solute and solvent must be present.
   if A & L, %#ok
      u(:,L) = u(:,L) + B(:,L).*R(:,V).*R(:,A).*z(:,A);
   end

   % If we allow mixing, both solute and solvent must be present.
   if A & V, %#ok
      u(:,V) = u(:,V) - B(:,V).*R(:,A).*z(:,A);
   end

   d      = ones(size(B, 1), 1) - R(:,L).*R(:,V);
   u      = bsxfun(@rdivide, u, d);
end

function rho = phaseDensity (B, R, rho_S, col)
   % rho = inv(B)*R'*rho_S, where rho_S is component surface density.
   A = col{1};  L = col{2};  V = col{3};
   if A, rho(:,A) = rho_S(A) ./ B(:,A); end
   if L, rho(:,L) = rho_S(L) ./ B(:,L); end
   if V, rho(:,V) = rho_S(V) ./ B(:,V); end

   % If we allow mixing, both solute and solvent must be present.
   if L & V, %#ok
      rho(:,L) = rho(:,L) + R(:,L)*rho_S(V)./B(:,L);
      rho(:,V) = rho(:,V) + R(:,V)*rho_S(L)./B(:,V);
   end

   % If we allow mixing, both solute and solvent must be present.
   if V & A, %#ok
      rho(:,A) = rho(:,A) + R(:,A)*rho_S(A)./B(:,A);
   end
end

function c   = phaseCompressibility (B, dB, R, dR, col)
   A = col{1};  L = col{2};  V = col{3};
   d        = ones(size(B, 1), 1) - R(:,L).*R(:,V);

   if A, c(:,A)   = -dB(:,A)./B(:,A); end
   if L, c(:,L)   = -dB(:,L)./B(:,L); end
   if V, c(:,V)   = -dB(:,V)./B(:,V); end

   % If we allow mixing, both solute and solvent must be present.
   if L & V, %#ok
      c(:,L) = c(:,L) + (B(:,V) - R(:,V).*B(:,L)).*dR(:,L)./(d.*B(:,L));
      c(:,V) = c(:,V) + (B(:,L) - R(:,L).*B(:,V)).*dR(:,V)./(d.*B(:,V));
   end

   % If we allow mixing, both solute and solvent must be present.
   if V & A, %#ok
      c(:,A) = c(:,A) + (B(:,V) - R(:,V).*B(:,L)).*dR(:,A)./(d.*B(:,A));
   end
end
