function mu = poise()
%Compute numerical value, in units of Pa*s, of one poise (P).
%
% SYNOPSIS:
%   mu = poise()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   mu - Numerical value, in units of Pascal seconds (Pa*s) of a viscosity
%        of one poise (P).

%{
#COPYRIGHT#
%}

% $Id: poise.m 1953 2009-03-31 10:54:12Z bska $

   mu = 0.1 * Pascal()*second();
end
