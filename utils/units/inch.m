function d = inch()
%Distance of one inch (in units of meters).
% 
% SYNOPSIS:
%   d = inch()
%
% PARMETERS:
%   None.
%
% RETURNS:
%   d - Distance of one inch (in units of meters).
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
#COPYRIGHT#
%}

% $Id: inch.m 1953 2009-03-31 10:54:12Z bska $

   d = 2.54 * centi()*meter();
end
