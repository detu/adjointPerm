function fluid = adjointFluidFields(fluid)
% add fields to fluid object that are needed in adjoint 

dLt  = @(sol)(fluid.dkr(sol)*(1./fluid.mu)');
d2Lt = @(sol)(fluid.d2kr(sol)*(1./fluid.mu)');

dLtInv  = @(sol)(-dLt(sol)./fluid.Lt(sol).^2);

d2LtInv = @(sol)(-(d2Lt(sol)./fluid.Lt(sol).^2 ...
                 -2*dLt(sol).^2./fluid.Lt(sol).^3));

D2fw = @(sol,kr, dkr, d2kr)( ...
            (1/fluid.mu(1))*(d2kr(:,1)./fluid.Lt(sol) - ...
            2*dkr(:,1).*dLt(sol)./fluid.Lt(sol).^2) ...
            + kr(:,1)./fluid.mu(1).*d2LtInv(sol));          


d2fw = @(sol)( D2fw(sol,fluid.kr(sol), fluid.dkr(sol), fluid.d2kr(sol)));

fluid.dLtInv  = dLtInv;
fluid.d2LtInv = d2LtInv;
fluid.d2fw    = d2fw;