function fluid = initSatnumFluid(T, varargin)
% Initialize incompressible two-phase fluid model with satnum and relperm
% evalution function.
%
% SYNOPSIS:
%   fluid = initSatnumFluid(T, varargin)
%
% PARAMETERS:
% PARAMETERS:
%   T       - Struct of (Eclipse) PVT tables.  It must contain
%             .swof.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - verbose -- Whether or not to emit informational messages
%                            while interpolating PVT tables.
%                            Logical.  Default value = FALSE.
%               - mu  -- Viscosities. Default value = [1 10]     [cP]
%               - rho -- Densities.   Default value = [1000 700] [kg/m^3]
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
%                    - relperm -- evaluates relative permeability curves
%                                 specified by SWOF.
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
%   initSatnumRelPerm, initFluid, initResSol, initWellSol, solveIncompFlow.

opt = struct(  'mu' , [1.0,     10.0], ...
               'rho', [1000.0, 700.0], ...
               'verbose', false);
opt = merge_options(opt, varargin{:});

assert(isfield(T, 'swof'))

relperm  = initSatnumRelPerm(T, 'verbose', opt.verbose);

   function kr   = return_kr(sol)
      kr = relperm(sol);
   end

   function dkr = return_dkr(sol)
      [trash, dkr] = relperm(sol);
   end

fluid = struct('N',   2,                                ...
               'kr',  @return_kr,                       ...
               'dkr', @return_dkr,                      ...
               'mu' , convertFrom(opt.mu, centi*poise), ...
               'rho', opt.rho);

fluid = initFluid(fluid);

% hack to replace relperm from call above:
fluid.relperm = relperm;

end
