function d = year()
%Give numerical value, in units of seconds, of one year.
%
% SYNOPSIS:
%   d = year()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   d - Numerical value, in units of s, of one year.
%
% SEE ALSO:
%   day, hour.

%{
#COPYRIGHT#
%}

% $Id: year.m 1953 2009-03-31 10:54:12Z bska $

   d = 365 * day();
end
