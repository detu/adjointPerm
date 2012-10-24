function [i, j, k] = getIndex(idx, ni, nj)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [i, j, k] = getIndex(idx, ni, nj)
%
% DESCRIPTION:
%   Find the i, j, k indices in a ni x nj x nk grid
%
% PARAMETERS:
%   idx         -   index from array of length (ni x nj x nk)
%
% RETURNS:
%   i, j, k     -   i, j, and k indices in grid with dimenions ni x nj x nk 
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
k = ceil(idx/(ni*nj));
j = ceil((idx-(k-1)*(ni*nj))/ni);           
i = ceil((idx-(k-1)*(ni*nj)-(j-1)*ni));