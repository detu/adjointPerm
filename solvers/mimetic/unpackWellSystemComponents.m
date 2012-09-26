function [BI, BIV, C, D, f, h] = unpackWellSystemComponents(W, mob)
%Extract hybrid linear system components from wells.
%
% SYNOPSIS:
%   [BIW, BIVW, CW, DW, fW, hW] = unpackWellSystemComponents(W, mob)
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
%   mob - Array, length n, of total mobility values.  One (positive) total
%         mobility value for each of the 'n' cells in the model.
%
% RETURNS:
%   BIW  - Cell array, one entry for each well, of "inverse mass
%          matrices" defined by well model indices.
%   BIVW - Cell array, one entry for each well, of matrix products
%          inv(B)*V, `V' denoting the flux in all perforations of a
%          single well.  Only relevant for compressible (i.e. black-oil)
%          flows (whenever field `BIV' is present in the well
%          structure), and all empty otherwise.
%   CW   - Cell array, one entry for each well, of well discrete
%          gradient oprators.
%   DW   - Cell array, one entry for each well, of well flux continuitiy
%          matrices.
%   fW   - Cell array, one entry for each well, of well target values.
%          Non-zero for pressure controlled (BHP) wells.  Zero otherwise.
%   hW   - Cell array, one entry for each well, of well target values.
%          Non-zero for rate controlled wells.  Empty otherwise.
%
% NOTE:
%   If the well structure is empty, i.e., if ISEMPTY(W), then all return
%   values are CELL(1) (a one-element cell array containing only an empty
%   numeric array).
%
% SEE ALSO:
%   solveIncompFlow, assembleWellSystem, unpackWellSystemComponentsMS.

%{
#COPYRIGHT#
%}

% $Id: unpackWellSystemComponents.m 1953 2009-03-31 10:54:12Z bska $

if isempty(W),
   BI  = { [] };
   BIV = { [] };
   C   = { [] };
   D   = { [] };
   f   = { [] };
   h   = { [] };
else
   S   = [ W.S   ];
   RHS = [ S.RHS ];

   numWells = numel(S);

   %---- System matrix components -------

   if isfield(S, 'BI'),
      BI = { S.BI };
   else
      BI = { [] };
   end

   if isfield(S, 'BIV'),
      BIV = { S.BIV };
   else
      BIV = cell(size(BI));
   end

   % If mobilities are provided by input:
   if ~isempty(mob) && ~isempty(BI{1}),
      for w = 1 : numWells,
         nc    = numel(W(w).cells);
         BI{w} = spdiags(mob(W(w).cells), 0, nc, nc) * BI{w};
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
end
