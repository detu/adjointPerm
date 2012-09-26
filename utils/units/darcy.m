function d = darcy()
%Compute numerical value, in units of m^2, of the Darcy constant.
%
% SYNOPSIS:
%   d = darcy()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   d - Numerical value, in units of m^2, of the Darcy constant.
%
% NOTE:
%   A porous medium with a permeability of 1 darcy permits a flow (flux) of
%   1 cm³/s of a fluid with viscosity 1 cP (1 mPa·s) under a pressure
%   gradient of 1 atm/cm acting across an area of 1 cm².
%
% SEE ALSO:
%   gravity.

%{
#COPYRIGHT#
%}

% $Id: darcy.m 1953 2009-03-31 10:54:12Z bska $

   mu  = centi()*poise(); % Fluid viscosity, 1cP = 1mPa·s
   cm  = centi()*meter(); % [m]
   sec = second();        % [s]

   p_grad = atm() / cm;   % Pressure gradient        [Pa/m]
   area   = cm ^ 2;       % Active area              [m^2]
   flow   = cm ^ 3 / sec; % Flow rate                [m^3/s]
   vel    = flow / area;  % Fluid velocity           [m/s]

   d = vel * mu / p_grad; % == 1e-7 [m^2] / 101325
                          % == 9.869232667160130e-13 [m^2]
end
