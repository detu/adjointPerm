function a = atm()
%Compute numerical value, in units of Pascal, of one atmosphere.
%
% SYNOPSIS:
%   a = atm()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   a - Numerical value, in units of Pascal (Pa), of a pressure of one
%       atmosphere.

%{
#COPYRIGHT#
%}

% $Id: atm.m 1953 2009-03-31 10:54:12Z bska $

   a = 101325 * Pascal();
end
