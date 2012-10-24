function [U] = fixStates(U,nstat)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [U] = fixStates(U,nstat)
%
% DESCRIPTION:
%   Ensure physical values for updates state variables
%
% PARAMETERS:
%   U           -   ensemble of model states
%   nstat       -   array containing number of elements for each state
%                   variables
%
% RETURNS:
%   U           -   physically valid state variable values
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
prme = U(1:nstat(1),:);
pore = U(sum(nstat(1:1))+1:sum(nstat(1:2)),:);
prfe = U(sum(nstat(1:2))+1:sum(nstat(1:3)),:);
sate = U(sum(nstat(1:3))+1:sum(nstat(1:4)),:);
if length(nstat) > 4
    stre = U(sum(nstat(1:4))+1:sum(nstat(1:5)),:);
else
    stre = [];
end

% porosity
it=find(pore<0.02); if ~isempty(it), pore(it)=0.02+(0.02-pore(it)).*rand(length(it),1); end
it=find(pore>0.40); if ~isempty(it), pore(it)=0.40-(pore(it)-0.40).*rand(length(it),1); end

U = [prme; pore; prfe; sate; stre];
