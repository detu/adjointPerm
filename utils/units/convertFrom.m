function q = convertFrom(q, unit)
%Convert physical quantity from given unit to equivalent SI.
%
% SYNOPSIS:
%   q = convertFrom(q, unit)
%
% PARAMETERS:
%   q    - Numerical array containing values of a physical quantity
%          measured in a given unit of measurement.
%
%   unit - The unit of measurement of the physical quantity 'q'.  Assumed
%          to be a combination of the known units in 'units'.
%
% RETURNS:
%   q    - Numerical array containing the numerical values resulting from
%          converting the input array from the unit of measurement given by
%          'unit' to the equivalent SI unit.
%
% EXAMPLES:
%   press = convertFrom(press, barsa())    % bar -> Pascal
%   rate  = convertFrom(rate, stb()/day()) % stb/day -> m^3/s
%   mu    = convertFrom(mu, Pa()*sec())    % Pa s -> Pa s (identity)
%
%   press = convertTo(convertFrom(press, atm()), psia()) % Atm -> psia
%
% NOTE:
%   It is the caller's responsibility to supply a 'unit' consistent with
%   the physical quantity 'q'.  Specifically, function 'convertFrom' does
%   no checking of the input parameters and will, if so instructed, convert
%   the numeric value of a pressure into a time value.  Caveat emptor.
%
% SEE ALSO:
%   units, convertTo.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

   q = q .* unit;
end
