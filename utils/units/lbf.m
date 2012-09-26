function F = lbf()
%Force excerted by a mass of one avoirdupois pound at Tellus equator.
% 
% SYNOPSIS:
%   F = lbf()
%
% PARMETERS:
%   None.
%
% RETURNS:
%   F - Force excerted by a mass of one avoirdupois pound (0.45359237 kg)
%       at Tellus equator measured in units of Newton (N).
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
#COPYRIGHT#
%}

% $Id: lbf.m 1953 2009-03-31 10:54:12Z bska $

   g = 9.80665;
   F = 0.45359237 * g;
end
