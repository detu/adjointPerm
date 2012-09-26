function fluid = initCoreyFluid(varargin)
%Initialize incompressible two-phase fluid model.
%
% SYNOPSIS:
%   fluid = initCoreyFluid
%   fluid = initCoreyFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - mu  -- Viscosities.  Default value = [1 10]     [cP]
%               - rho -- Densities.    Default value = [1000 700] [kg/m^3]
%               - n   -- Exponents.    Default value = [2, 2]
%               - sr  -- Residual sat. Default value = [0.2, 0.2]
%               - kwm -- kwm=kr(sr).   Default value = [1, 1]
%
% RETURNS:
%   fluid - Fluid data structure representing the current state of the
%           fluids within the reservoir model. Contains scalar fields and
%           function handles that take a structure sol containing a field s
%           (saturation) as argument. sol is normally the reservoir
%           solution structure.
%               -- Scalar fields:
%                    - N      -- Number of phases = 2.
%                    - mu     -- Viscosities.
%                    - rho    -- Densities.
%               -- Function handles:
%                    - kr     -- Relative permeabilities.
%                    - dkr    -- Derivative of relative permeabilities.
%                    - mob    -- Mobilities.
%                    - Lt     -- Total mobility.
%                    - omega  -- Accumulated phase densities weighted by
%                                fractional flow functions.
%
% NOTE:
%   The fluid model represented by the return argument is the two-phase
%   incompressible counterpart to the fluid model of the Black Oil 'pvt'
%   function. For generality, the argument to all member functions is
%   assumed to be a structure that contains (at least) a saturation field,
%   's'.  Therefore, to plot the relative permeability, use statements akin
%   to
%
%      s = linspace(0, 1, 20).'; kr = fluid.kr(struct('s', s));
%      plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   initSimpleFluid, initFluid, initResSol, initWellSol, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2009-11-24 18:41:11 +0100 (Tue, 24 Nov 2009) $
% $Revision: 3222 $

   opt = struct('mu' ,    [1.0,     10.0], ...
                'rho',    [1000.0, 700.0], ...
                'n'  ,    [2,          2], ...
                'sr' ,    [0.2,      0.2], ...
                'kwm',    [1,          1]);

   opt = merge_options(opt, varargin{:});

   kr    = @(sol) relperm  (sol, opt);
   dkr   = @(sol) drelperm (sol, opt);
   d2kr  = @(sol) d2relperm(sol, opt);

   
   fluid = struct('N'   , 2,                                ...
                  'kr'  , kr,                               ...
                  'dkr' , dkr,                              ...
                  'd2kr', d2kr,                             ...
                  'mu'  , convertFrom(opt.mu, centi*poise), ...
                  'rho' , opt.rho,                          ...
                  'sr'  , opt.sr);
   fluid = initFluid(fluid);

end

%--------------------------------------------------------------------------
% Helpers (i.e., return values) follow.
%--------------------------------------------------------------------------

function q = relperm(sol, opt)
   [s1, s2] = modified_saturations(sol, opt);
   q        = [opt.kwm(1) * s1.^opt.n(1), opt.kwm(2)*s2.^opt.n(2)];
end

%--------------------------------------------------------------------------

function q = drelperm(sol, opt)
   [s1, s2, den] = modified_saturations(sol, opt);
   n = opt.n;
   q = [ opt.kwm(1) * n(1).*s1.^(n(1)-1), ...
        -opt.kwm(2) * n(2).*s2.^(n(2)-1)] ./ den;
end

%--------------------------------------------------------------------------

function q = d2relperm(sol, opt)
   [s1, s2, den] = modified_saturations(sol, opt);
   n = opt.n;
   q = [ opt.kwm(1) * n(1)*(n(1)-1).*s1.^(n(1)-2), ...
         opt.kwm(2) * n(2)*(n(2)-1).*s2.^(n(2)-2)] ./ den.^2;
end


%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(sol, opt)
   den      = 1 - sum(opt.sr);
   s1       = (    sol.s(:,1) - opt.sr(1)) ./ den;
   s1(s1<0) = 0;
   s2       = (1 - sol.s(:,1) - opt.sr(2)) ./ den;
   s2(s2<0) = 0;
end
