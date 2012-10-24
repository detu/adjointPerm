function [grid,well,fluid,rock,p0,s0,phi,K,G] = getModel(workdir,modelname,ne,j,G)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [grid,well,fluid,rock,p0,s0,phi,K,G] = getModel(workdir,modelname,ne,j,G)
%
% DESCRIPTION:
%   Obtain a reservoir model description including grid structure, fluid
%   and rock properties, well descriptions, and the initial state, as well
%   as ensemble of permeability, porosity and structural variables (if
%   structural updating is desired).
%
% PARAMETERS:
%   workdir     -   working folder
%   modelname   -   name of the model (there has to be a .m file with the
%                   same name)
%   ne, j       -   ensemble size and current member
%   G           -   structural parameter array
%
% RETURNS:
%   grid        -   structural grid
%   well        -   description of well in the model
%   fluid       -   fluid properties
%   rock        -   permeability and porosity of the truth model
%   p0, s0      -   initial pressure and saturation
%   phi, K, G   -   ensembles of porosity, permeability and structure
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
if exist([workdir modelname '.mat'],'file') && isempty(G)
    load([modelname '.mat'],'grid','well','fluid','rock','p0','s0','phi','K');
else
    [grid,well,fluid,rock,p0,s0,phi,K,G] = makeModel(modelname,ne,j,G);
    save([modelname '.mat'],'grid','well','fluid','rock','p0','s0','phi','K','G');
end

return