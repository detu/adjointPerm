function [B, C, D, f, h] = unpackMixedWellSystemComponents(W, mob)
% unpackWellSystemComponents -- Extract hybrid linear system components
%                               from wells.
%
% SYNOPSIS:
%   [BW, CW, DW, fW, hW] = unpackMixedWellSystemComponents(W, mob)
%
% DESCRIPTION:
%   This function recovers contributions to the final, hybrid, system of
%   linear equations from which velocities and cell as well contact
%   (face) pressures are determined.
%
% PARAMETERS:
%   W   - Well structure as defined by functions addWell,
%         assembleWellSystem, and (e.g.) updateBlackOilWellSystem.
%
%   mob - List of propertyname/propertyvalue pairs.
%             PropertyNames     PropertyValues     Default
%             'MobilityValues'  vector             ones
%
% RETURNS:
%   BW   - Cell array, one entry for each well, of "mass matrices" defined
%          by well model productivity indices.
%   CW   - Cell array, one entry for each well, of well discrete
%          gradient oprators.
%   DW   - Cell array, one entry for each well, of well flux continuitiy
%          matrices.
%   fW   - Cell array, one entry for each well, of well target values.
%          Non-zero for pressure controlled (BHP) wells.  Zero otherwise.
%   hW   - Cell array, one entry for each well, of well target values.
%          Non-zero for rate controlled wells.  Empty otherwise.

% $Id: unpackMixedWellSystemComponents.m 1892 2009-03-25 16:38:40Z bska $

S   = [ W.S   ];
RHS = [ S.RHS ];

numWells = length(S);

%---- System matrix components -------

B = { S.B };

% If mobilities are provided by input:
if ~isempty(mob),
   for w = 1 : numWells,
      nc   = numel(W(w).cells);
      B{w} = spdiags(1 ./ mob(W(w).cells), 0, nc, nc) * B{w};
   end
end

C = { S.C };
D = { S.D };

%---- Right hand side components -----

if isfield(RHS, 'f_grav'),
   f = cellfun(@plus, { RHS.f }, { RHS.f_grav }, ...
               'UniformOutput', false);
else
   f = { RHS.f };
end
h  = { RHS.h };
