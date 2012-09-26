function [ff, gg, hh, gp, dF, dC] = computePressureRHS(g, omega, bc, src)
%Compute right-hand side contributions to pressure linear system.
%
% SYNOPSIS:
%   [f, g, h, grav, dF, dC] = computePressureRHS(G, omega, bc, src)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   omega - Accumulated phase densities \rho_i weighted by fractional flow
%           functions f_i -- i.e., omega = \sum_i \rho_i f_i.  One scalar
%           value for each cell in the discretised reservoir model, G.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Pressure (f), source/sink (g), and flux (h) external
%             conditions.  In a problem without effects of gravity, these
%             values may be passed directly on to linear system solvers
%             such as 'schurComplementSymm'.
%
%   grav    - Pressure contributions from gravity.  One scalar value for
%             each half-face in the model (size(G.cellFaces,1)).
%
%   dF, dC  - Packed Dirichlet/pressure condition structure.  The faces in
%             'dF' and the values in 'dC'.  May be used to eliminate known
%             face pressures from the linear system before calling a system
%             solver (e.g., 'schurComplementSymm').
%
% SEE ALSO:
%   addBC, addSource, pvt, computeMimeticIP, schurComplementSymm.

%{
#COPYRIGHT#
%}

% $Date: 2009-07-01 12:26:03 +0200 (on, 01 jul 2009) $
% $Revision: 2421 $

   gp = grav_pressure(g, omega);
   ff = zeros(size(gp));
   gg = zeros([g.cells.num, 1]);
   hh = zeros([g.faces.num, 1]);

   if ~isempty(src),
      assert (max([src.cell]) <= g.cells.num);

      ss = accumarray(src.cell, src.rate)    ;
      ii = accumarray(src.cell,    1    ) > 0;
      gg(ii) = gg(ii) + ss(ii);
   end

   dF = false([g.faces.num, 1]);
   dC = [];

   if ~isempty(bc),
      assert (max([bc.face]) <= g.faces.num);
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      % Pressure (Dirichlet) boundary conditions.
      %  1) Extract the faces marked as defining pressure conditions.
      %     Define a local numbering (map) of the face indices to the
      %     pressure condition values.
      %
      is_press = strcmpi('pressure', bc.type);
      face     = bc.face (is_press);
      dC       = bc.value(is_press);
      map      = sparse(double(face), 1, 1 : numel(face));

      %  2) For purpose of (mimetic) pressure solvers, mark the 'face's as
      %     having pressure boundary conditions.  This information will be
      %     used to eliminate known pressures from the resulting system of
      %     linear equations.  See (e.g.) 'solveIncompFlow' for details.
      %
      dF(face) = true;

      %  3) Enter Dirichlet conditions into system right hand side.
      %     Relies implictly on boundary faces being mentioned exactly once
      %     in g.cellFaces(:,1).
      %
      i     =   dF(    g.cellFaces(:,1) );
      ff(i) = - dC(map(g.cellFaces(i,1)));

      %  4) Reorder Dirichlet conditions according to SORT(face).  This
      %     allows the caller to issue statements such as 'X(dF) = dC' even
      %     when ISLOGICAL(dF).
      %
      dC = dC(map(dF));

      % Flux (Neumann) boundary conditions.
      % Note negative sign due to bc.value representing INJECTION flux.
      %
      is_flux              = strcmpi('flux', bc.type);
      hh(bc.face(is_flux)) = - bc.value(is_flux);
   end

   assert (~any(dC < 0));  % Pressure conditions should always be non-neg.
end

%--------------------------------------------------------------------------

function ff = grav_pressure(g, omega)
   g_vec = gravity();   % Must be a 1-by-3 row vector for subsequent code.
   dim   = size(g.nodes.coords,2);

   assert (1 < dim && dim < 4);
   assert (all(size(g_vec) == [1,3]));

   if norm(g_vec(1 : dim)) > 0,
      cellno = rldecode(1 : g.cells.num, double(g.cells.numFaces), 2) .';
      cvec   = g.faces.centroids(g.cellFaces(:,1), :) - ...
               g.cells.centroids(cellno          , :);
      ff     = omega(cellno) .* (cvec * g_vec(1:dim).');
   else
      ff     = zeros([size(g.cellFaces,1), 1]);
   end
end
