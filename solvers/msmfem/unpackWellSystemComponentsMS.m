function [B, C, D, f, h, Psi, Phi] = ...
      unpackWellSystemComponentsMS(W, G, p, mob, omega)
%Extract coarse linear system components from wells.
%
% SYNOPSIS:
%   [B, C, D, f, h, ...
%    Psi, Phi] = unpackWellSystemComponentsMS(W, G, p, mob, rho)
%
% DESCRIPTION:
%   This function recovers contributions to the final, hybrid, system of
%   linear equations from which velocities and cell as well as contact
%   (face) pressures are determined.
%
% PARAMETERS:
%   W   - Well structure as defined by functions addWell and
%         generateCoarseWellSystem.
%
%   G   - Underlying fine grid model of reservoir.
%
%   p   - Cell-to-block partition vector.
%
%   mob - Array, length n, of total mobility values.  One (positive) total
%         mobility value for each of the 'n' cells in the model.
%
%   rho - Array, length nf, of fluid densities.  One (positive) density
%         value for each fluid in the simulation model.  Only used in the
%         presence of gravity effects to add hydrostatic pressure
%         adjustments within the wells.
%
% RETURNS:
%   B   - Cell array, one entry for each well, of "mass matrices"
%         defined by well model productivity indices.
%
%   C   - Cell array, one entry for each well, of well discrete gradient
%         operators.
%
%   D   - Cell array, one entry for each well, of well flux continuitiy
%         matrices.
%
%   f   - Cell array, one entry for each well, of well target values.
%         Non-zero for pressure controlled (BHP) wells.  Zero otherwise.
%
%   h   - Cell array, one entry for each well, of well target values.
%         Non-zero for rate controlled wells.  Empty otherwise.
%
%   Psi - Well flux basis functions represented as a single sparse matrix
%         of size nfhf-by-(sum_w (numb of coarse blocks intersected by w)).
%         'nfhf' is the number of fine-scale half-faces.
%
%   Phi - Well pressure basis functions represented as a single sparse
%         matrix of size G.cells.num-by-(sum_w (above)).
%
% NOTE:
%   If the well structure is empty, i.e., if ISEMPTY(W), then all return
%   values, save 'Psi' and 'Phi', are CELL(1) (a one-element cell array
%   containing only an empty numeric array).  In this case 'Psi' and 'Phi'
%   are empty sparse matrices of the correct size (correct number of rows
%   and no columns).
%
% SEE ALSO:
%   solveIncompFlowMS, gravity, generateCoarseWellSystem, evalWellBasis,
%   unpackWellSystemComponents.

%{
#COPYRIGHT#
%}

% $Id: unpackWellSystemComponentsMS.m 2338 2009-06-05 17:19:30Z bska $

   if isempty(W),
      [B, C, D, f, h] = deal({ [] });
      Psi = sparse(size(G.cellFaces,1), 0);
      Phi = sparse(G.cells.num        , 0);
   else
      CS  = [ W.CS   ];
      RHS = [ CS.RHS ];

      numWells = length(CS);

      %---- System matrix components -------

      use_grav = norm(gravity()) > 0;
      if use_grav,
         % Expand well pressure rhs to each perforation, adjust for gravity
         % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
         nperf  = cellfun(@numel, { W.cells });
         i      = 0;

         dp     = norm(gravity()) * vertcat(W.dZ) .* ...
                  omega(vertcat(W.cells));
         f_grav = cell([1, numWells]);
      end

      B  = cell([numWells, 1]);
      for w = 1 : numWells,
         nc    = numel(W(w).cells);
         rates = W(w).CS.rates;

         mult = sparse(1:nc, 1:nc, mob(W(w).cells) .* W(w).WI, nc, nc);
         B{w} = rates.' * (mult \ rates);

         if use_grav,
            f_grav{w} = rates.' * dp(i + 1 : i + nperf(w));
            i         = i + nperf(w);
         end
      end

      % Form flux basis function matrix by extracting contributions for all
      % wells.  The SPARSE input triplets are stored directly in the first
      % three entries of the block '.basis' field for each well (details in
      % function 'evalWellBasis'), so simply extract these values.  We
      % build a (large) cell array of these cell arrays by concatenating
      % this field for all wells.  Note that the 'j' vectors contain only a
      % single 1-entry in j(1), so we form the required input vector 'j' by
      % computing the cumulative sum of these values (i.e., put the values
      % into consecutive columns of 'Psi').
      [i,j,v] = cellfun(@(x) x{1:3}, vertcat(CS.basis), ...
                        'UniformOutput', false);
      Psi     = sparse(       vertcat(i{:}) , ...
                       cumsum(vertcat(j{:})), ...
                              vertcat(v{:}) , ...
                       size(G.cellFaces,1), numel(j));

      % Form pressure basis function matrix by extracting contributions for
      % all wells.  The SPARSE input triplets are stored directly in the
      % first three entries of the block '.basisP' field for each well, so
      % simply extract these values.
      [i,j,v] = cellfun(@(x) x{1:3}, vertcat(CS.basisP), ...
                        'UniformOutput', false);
      Phi     = sparse(       vertcat(i{:}) , ...
                       cumsum(vertcat(j{:})), ...
                              vertcat(v{:}) , G.cells.num, numel(j));

      % Update the pressure basis function for effects of (modified)
      % average mobility.
      %
      % 1) Compute current volume-weighted average total mobility in each
      %    coarse block.
      lam = accumarray(p, mob .* G.cells.volumes) ./ ...
            accumarray(p,        G.cells.volumes);

      % 2) Extract average block mobilities by which the individual basis
      %    functions were generated.
      avgmob = reshape(cellfun(@(x) x{4}, vertcat(CS.basisP)), [], 1);

      % The well pressure basis functions are entered into the '.basisP'
      % cell arrays in order of *increasing* coarse block numbers (function
      % 'evalWellBasis' calls UNIQUE before entering the evaluation loop on
      % line 82).  Therefore, the correct average mobility update factor is
      % obtained simply by ordering the unique coarse blocks intersected by
      % a well in order of increasing block numbers.
      %
      pw  = cellfun(@(x) unique(p(x)), { W.cells }, ...
                    'UniformOutput', false);

      % Perform actual pressure basis function update.
      nb  = size(Phi,2);
      Phi = Phi * spdiags(avgmob ./ lam(vertcat(pw{:})), 0, nb, nb);

      % Extract remaining coarse (well) system components.
      C = { CS.C  };
      D = { CS.D  };
      f = { RHS.f };
      h = { RHS.h };

      if use_grav,
         f = cellfun(@minus, f, f_grav, 'UniformOutput', false);
      end
   end
end
