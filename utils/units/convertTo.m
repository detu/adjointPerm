function q = convertTo(q, unit)
%Convert physical quantity from SI to equivalent given unit.
%
% SYNOPSIS:
%   q = convertTo(q, unit)
%
% PARAMETERS:
%   q    - Numerical array containing values of a physical quantity
%          measured in a strictly SI unit.
%
%   unit - The unit of measurement to which 'q' should be converted.
%          Assumed to be a combination of the known units in 'units'.
%
% RETURNS:
%   q    - Numerical array containing the numerical values resulting from
%          converting the input array from the SI unit to the unit of
%          measurement given by 'unit'.
%
% EXAMPLES:
%   press = convertTo(press, barsa())    % Pascal -> bar
%   rate  = convertTo(rate, stb()/day()) % m^3/s -> stb/day
%   mu    = convertTo(mu, Pa()*sec())    % Pa s -> Pa s (identity)
%
%   press = convertTo(convertFrom(press, atm()), psia()) % Atm -> psia
%
% NOTE:
%   It is the caller's responsibility to supply a 'unit' consistent with
%   the physical quantity 'q'.  Specifically, function 'convertTo' does
%   no checking of the input parameters and will, if so instructed,
%   convert the numeric value of a pressure into a time value.  Caveat
%   emptor.
%
% SEE ALSO:
%   units, convertFrom.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

   q = q ./ unit;
end
