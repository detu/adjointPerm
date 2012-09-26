function t = hour()
%Time span of one hour (in units of seconds).
% 
% SYNOPSIS:
%   t = hour()
%
% PARMETERS:
%   None.
%
% RETURNS:
%   t - Time span of one hour measured in units of seconds.
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
#COPYRIGHT#
%}

% $Id: hour.m 1953 2009-03-31 10:54:12Z bska $

   t = 60 * minute();
end
