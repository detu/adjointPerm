function fluid = initFluid(fluid)
%Complete initialization of incompressible two-phase fluid model.
%
% SYNOPSIS:
%   fluid = initFluid(fluid)
%
% PARAMETERS:
%   fluid - fluid object to be processed.
%
% RETURNS:
%   fluid - Fluid data structure representing the current state of the
%           fluids within the reservoir model.
%
% NOTE:
%   The fluid model represented by the return argument is the two-phase
%   incompressible counterpart to the fluid model of the Black Oil 'pvt'
%   function.
%
% SEE ALSO:
%   initSimpleFluid, initResSol, initWellSol, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
% $Revision: 2926 $

   mob = @(sol)   bsxfun(@rdivide, fluid.kr(sol), fluid.mu);
   w1  = @(mob,e) (bsxfun(@times, fluid.rho, mob) * e) ./ (mob * e);
   Lt  = @(mob)      mob * ones([size(mob,2), 1]) ;
   w   = @(mob)   w1(mob , ones([size(mob,2), 1]));
   
   fluid.mob   = mob;
   fluid.Lt    = @(sol) Lt(mob(sol));
   fluid.omega = @(sol) w (mob(sol));
   
   % add fractional flow function and derivative
   if fluid.N == 2
      c = 1;      
      if isfield(fluid, 'sr')
         c = 1/(1-sum(fluid.sr));
      end           
      fw  = @(kr) (kr(:,1)./ ...
                  (kr(:,1) + (fluid.mu(1)./fluid.mu(2))*kr(:,2)));
      dfw = @(kr, dkr, Lt) c*((dkr(:,1)./fluid.mu(1))./Lt - ...
                          (kr(:,1)./fluid.mu(1)).* ...      
                          sum(bsxfun(@rdivide, dkr, fluid.mu),2)./(Lt.^2));
      fluid.fw  = @(sol) fw(fluid.kr(sol));
      fluid.dfw = @(sol) dfw(fluid.kr(sol), fluid.dkr(sol), Lt(mob(sol)));
   end   

   function [mu, rho, c, u, R] = pvt (p, x)
      [c, u, R] = deal(nan);
      assert (numel(p) == size(x, 1));

      if     isa(fluid.mu, 'function_handle'),
         mu = fluid.mu(p);
      elseif isnumeric(fluid.mu),
         mu = repmat(fluid.mu, [numel(p), 1]);
      end

      if     isa(fluid.rho, 'function_handle'),
         rho = fluid.rho(p);
      elseif isnumeric(fluid.rho),
         rho = repmat(fluid.rho, [numel(p), 1]);
      end
   end
   fluid.pvt = @pvt;

   function [kr, dkr] = relperm (s)
      if     isa(fluid.kr, 'function_handle'),
         kr = fluid.kr(s);
      elseif isnumeric(fluid.kr),
         kr = repmat(fluid.kr, [size(s, 1), 1]);
      end

      if     isa(fluid.dkr, 'function_handle'),
         dkr = fluid.dkr(s);
      elseif isnumeric(fluid.dkr),
         dkr = repmat(fluid.dkr, [size(s, 1), 1]);
      end
   end
   fluid.relperm = @relperm;
end
