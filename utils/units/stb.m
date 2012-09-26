function bbl = stb()
%Compute numerical value, in units of m^3, of one standard barrel.
%
% SYNOPSIS:
%   bbl = stb()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   bbl - Numerical value, in units of m^3, of a volume of one barrel.

%{
#COPYRIGHT#
%}

% $Id: stb.m 1953 2009-03-31 10:54:12Z bska $

   bbl = 0.159 * meter()^3;
end
