function fluid = initTwoPhaseIncompressibleFluid(varargin)
% initTwoPhaseIncompressibleFluid -- Initialize incompressible two-phase fluid model
%
% SYNOPSIS:
%   fluid = initTwoPhaseIncompressibleFluid
%   fluid = initTwoPhaseIncompressibleFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - muw  -- Water viscosity.            Default value =    0.3
%               - rhow -- Water density.              Default value = 1000.0
%               - muo  -- Oil viscosity.              Default value =    2.85
%               - rhoo -- Oil density.                Default value =  700.0
%               - swc  -- Connate water saturation    Default value =    0
%               - sor  -- Irreducible oil saturation  Defailt value =    0
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
%   initSimpleFluid, pvt, initResSol, initWellSol, solveWellSystem.


opt = struct('muw' ,    0.3 , ...
             'muo' ,    2.85, ...
             'rhow', 1000.0 , ...
             'rhoo',  700.0 , ...
             'swc' ,    0   , ...
             'sor' ,    0   );
opt = merge_options(opt, varargin{:});

muw   = opt.muw;                             assert (muw > 0);
muo   = opt.muo;                             assert (muo > 0);
rhow  = opt.rhow;
rhoo  = opt.rhoo;
swc   = opt.swc;
sor   = opt.sor;

ss    = @(s) ( s - swc )/( 1 - sor - swc );         % scaled saturatoin
Lw    = @(s) ( ss(s).^2 )/muw;                      % water mobility
Lo    = @(s) ( (1 - ss(s)).^2 )/muo;                % oil mobility
lt    = @(s) Lw(s) + Lo(s);
Lt    = @(sol) lt(sol.sw);                          % total mobility    
fw    = @(s) Lw(s)./lt(s);                          % fractional flow
omega = @(sol) (rhow*Lw(sol.sw) + rhoo*Lo(sol.sw))./Lt(sol);

% Derivatives:
mr    = muw/muo;
c     = 1/( 1 - sor - swc );
% fractional flow function
Dfw   = @(s) c*( (2*mr) * ss(s) .*(1-ss(s)) ) ./ ...
               ( ( ss(s).^2 + mr*(1-ss(s)).^2 ).^2 );


% Dfw   = @(s) ( 2*ss(s)*c ./ ( (ss(s).^2) + mr*((1-ss(s)).^2) ) ) - ...
%              ( (( ((2/muw)*ss(s)*c ) - (2/muo)*((1-ss(s))*c) ) .*(ss(s).^2) ) ./ ( muw*(ss(s).^2) + (muw*mr*((1-ss(s)).^2)) ).^2 );

% second order derivative
% D2fw = @(s) ( (4*c*mr* ( ( 2*ss(s)*c ) - 2*mr*c*(1 - ss(s)) ) .* (1 - ss(s)) .* ss(s) ) ./ (( ( ss(s).^2 ) + mr*((1-ss(s)).^2)  ).^3) ) + ...
%            ( 2*c*mr* (1-ss(s))* c ./ ( (ss(s).^2) + mr*((1-ss(s)).^2) ).^2 ) - ...
%            ( (2*c*mr*ss(s)*c) ./ ( ((ss(s).^2) + (mr*((1-ss(s)).^2) )).^2 ) );


D2fw = @(s) c*( (2*c*mr* ((1-ss(s)) - ss(s))  .* ( ss(s).^2 + mr*(1-ss(s)).^2) ) - ( 4*mr*ss(s).* (1-ss(s)) .*(2*ss(s)*c - 2*mr*(1-ss(s))*c )) ) ... 
            ./ ( ( ss(s).^2 + mr*(1-ss(s)).^2 ).^3 );


%inverse of mobility
DLtInv= @(s) c*( (-2/muw)*ss(s) + (2/muo)*(1-ss(s)) ) ./ ...
            ( ( (1/muw)*ss(s).^2 + (1/muo)*(1-ss(s)).^2 ).^2 );

% D2LtInv= @(s) ( -c^2*( (2/muw)  + (2/muo)  ) ./ ( ( ( (1/muw)* (ss(s).^2) ) + ( (1/muo)*((1-ss(s)).^2) ) ).^2) ) - ...
%               ( 8*c^2* ( ( ((1/muw)* ss(s)) - ( (1/muo)* (1-ss(s)) ) ) .* ( (1/muo)* (1-ss(s)) - ((1/muw)* ss(s)) ) ) ./ (( ((1/muw)* ss(s).^2) + (1/muo)*((1-ss(s)).^2) ).^3) );


D2LtInv= @(s) c*( ( -2*c*( (1/muw)+(1/muo) ) .* ( (1/muw)*ss(s).^2 + (1/muo)*(1-ss(s)).^2 ) ) - ( c*( ((-2/muw)*ss(s)) + ((2/muo)*(1-ss(s))) ) .* (4*c*( ((1/muw)*ss(s)) - ((1/muo)*(1-ss(s))) ) ) ))...
              ./ ( ( (1/muw)*ss(s).^2 + (1/muo)*(1-ss(s)).^2 ).^3 );


fluid = struct ('muw',   muw,  ...
                'muo',   muo,  ...
                'rhow',  rhow, ...
                'rhoo',  rhoo, ...
                'swc',   swc,  ...
                'sor',   sor,  ...
                'Lw',    Lw,   ...
                'Lo',    Lo,   ...
                'lt',    lt,   ...
                'Lt',    Lt,   ...
                'fw',    fw,   ...
                'omega', omega,...
                'Dfw',   Dfw,  ...
                'D2fw',  D2fw, ...
                'DLtInv', DLtInv, ...
                'D2LtInv', D2LtInv);
