function fluid = initSimpleFluid(varargin)
%Initialize incompressible two-phase fluid model.
%
% SYNOPSIS:
%   fluid = initSimpleFluid
%   fluid = initSimpleFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - mu  -- Viscosities. Default value = [1 10]     [cP]
%               - rho -- Densities.   Default value = [1000 700] [kg/m^3]
%               - n   -- Exponents.   Default value = [2, 2]
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
%   function.  For generality, the argument to all member functions is
%   assumed to be a structure that contains (at least) a saturation field,
%   's'.  Therefore, to plot the relative permeability, use statements akin
%   to
%
%      s = linspace(0, 1, 20).'; kr = fluid.kr(struct('s', s));
%      plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   initFluid, initResSol, initWellSol, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
% $Revision: 2926 $

opt = struct('mu' , [1.0,     10.0], ...
             'rho', [1000.0, 700.0], ...
             'n'  , [2,          2]);

opt = merge_options(opt, varargin{:});

n     = opt.n;
kr    = @(sol) [sol.s(:,1).^n(1), (1-sol.s(:,1)).^n(2)];
dkr   = @(sol) [ n(1) .*  sol.s(:,1)   .^(n(1)-1), ...
                -n(2) .* (1-sol.s(:,1)).^(n(2)-1)];

fluid = struct('N',   2,                                ...
               'kr',  kr,                               ...
               'dkr', dkr,                              ...
               'mu' , convertFrom(opt.mu, centi*poise), ...
               'rho', opt.rho);

fluid = initFluid(fluid);
