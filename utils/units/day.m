function d = day()
%Give numerical value, in units of seconds, of one day.
%
% SYNOPSIS:
%   d = day()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   d - Numerical value, in units of s, of one day.
%
% SEE ALSO:
%   darcy.

%{
#COPYRIGHT#
%}

% $Id: day.m 1953 2009-03-31 10:54:12Z bska $

   d = 24 * hour();
end
