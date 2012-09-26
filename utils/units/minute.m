function t = minute()
%Time span of one minute (in units of seconds).
% 
% SYNOPSIS:
%   t = minute()
%
% PARMETERS:
%   None.
%
% RETURNS:
%   t - Time span of one minute measured in units of seconds.
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
#COPYRIGHT#
%}

% $Id: minute.m 1953 2009-03-31 10:54:12Z bska $

   t = 60 * second();
end
