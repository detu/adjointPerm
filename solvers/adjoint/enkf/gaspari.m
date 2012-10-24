function [c0] = gaspari (z,c)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [c0] = gaspari (z,c)
%
% DESCRIPTION:
%   Evaluate the Gaspari-Cohn correlation function with local support
%
% PARAMETERS:
%   z           -   distance
%   c           -   correlation scale
%
% RETURNS:
%   c0          -   localization factor
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
if z < c
  c0=-0.25*(z/c)^5+0.5*(z/c)^4+0.625*(z/c)^3-(5.0/3.0)*(z/c)^2+1;
elseif z < 2*c
  c0=(1.0/12.0)*(z/c)^5-0.5*(z/c)^4+0.625*(z/c)^3+(5.0/3.0)*(z/c)^2-5*(z/c)+4-(2.0/3.0)*(c/z);
else
  c0=0;
end

return