function p = psia()
%Compute numerical value, in units of Pascal, of one Psi.
% 
% SYNOPSIS:
%   p = psia()
%
% PARMETERS:
%   None.
%
% RETURNS:
%   p - Numerical value, in units of Pascal (Pa), of a pressure of one Psi.
%
% NOTE:
%   The primary purpose of this utility function is to make examples
%   easier to read.

%{
#COPYRIGHT#
%}

% $Id: psia.m 1953 2009-03-31 10:54:12Z bska $

   p = lbf() / inch()^2;
end
