function [G,W,fluid,rock,p0,s0,phi,K,S] = makeModel(modelname,ne,j,S)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [G,W,fluid,rock,p0,s0,phi,K,S] = makeModel(modelname,ne,j,S)
%
% DESCRIPTION:
%   Obtain a reservoir model description including grid structure, fluid
%   and rock properties, well descriptions, and the initial state, as well
%   as ensemble of permeability, porosity and structural variables (if
%   structural updating is desired).
%
% PARAMETERS:
%   modelname   -   name of the model (there has to be a .m file with the
%                   same name)
%   ne, j       -   ensemble size and current member
%   S           -   structural parameter array
%
% RETURNS:
%   G           -   structural grid
%   W           -   description of well in the model
%   fluid       -   fluid properties
%   rock        -   permeability and porosity of the truth model
%   p0, s0      -   initial pressure and saturation
%   phi, K, S   -   ensembles of porosity, permeability and structure
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
if ~exist([modelname '.m'],'file')
    
      error('*** no valid model specified ***');
      
else
    
      eval(modelname)
      
end
