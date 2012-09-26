function [xr, xw] = incompMPFA(xr, xw, g, T, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   [xr, xw] = incompTPFA(xr, xw, G, T, fluid)
%   [xr, xw] = incompTPFA(xr, xw, G, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a multi-point flux approximation (MPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%   xr, xw - Reservoir and well solution structures either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompTPFA' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput - Whether or not to return the final system matrix 'A' to
%                  the caller of function 'incompTPFA'.
%                  Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%          - cellFlux     -- Outgoing flux across each local interface for
%                            all cells in the model.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'incompTPFA:DrivingForce:Missing'
%
% EXAMPLE:
%    G   = computeGeometry(cartGrid([3,3,5]));
%    f   = initSingleFluid();
%    rock.perm = rand(G.cells.num, 1)*darcy()/100;
%    bc  = pside([], G, 'LEFT', 1:G.cartDims(2), 1:G.cartDims(3),2);
%    src = addSource([], 1, 1);
%    W   = verticalWell([], G, rock, 1, G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'rate', 'Val', 1/day(), ...
%                       'InnerProduct', 'ip_tpf');
%    W   = verticalWell(W, G, rock, G.cartDims(1),   G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'bhp', ...
%                       'Val',  1*barsa(), 'InnerProduct', 'ip_tpf');
%    T   = computeMultiPointTrans(G, rock);
%    xr  = initResSol (G, 10);
%    xw  = initWellSol(G, 10);
%    [xr,xw] = incompMPFA(xr, xw, G, T, f, 'bc',bc,'src',src,'wells',W,...
%                         'MatrixOutput',true);
%
%    plotCellData(G, xr.cellPressure);
%
% SEE ALSO:
%   computeMultiPointTrans, addBC, addSource, addWell, initSingleFluid,
%   initResSol, initWellSol.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-26 15:26:29 +0200 (fr, 26 jun 2009) $
% $Revision: 2410 $


   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false);
   opt = merge_options(opt, varargin{:});

   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), gravity() == 0]),
      warning(id('DrivingForce:Missing'),                       ...
              ['No external driving forces present in model--', ...
               'state remains unchanged.\n']);
   end
   cellNo = rldecode(1:g.cells.num, double(g.cells.numFaces), 2) .';
   C      = sparse(1:size(g.cellFaces, 1), cellNo, 1);
   
   ind    = [(1:g.cells.num)'; max(g.faces.neighbors')'];
   totmob = fluid.Lt(xr);
   totmob = spdiags(totmob(ind), 0, numel(ind),numel(ind));
   T      = T*totmob;

   % Boundary conditions and source terms.
   % Note: Function 'computeRHS' is a modified version of function
   % 'computePressureRHS' which was originally written to support the
   % hybrid mimetic method.
   [ff, gg, hh, grav, dF, dC] = computeRHS(g, fluid.omega(xr), ...
                                           opt.bc, opt.src);
   b   = any(g.faces.neighbors==0, 2);
   I1  = [(1:g.cells.num)';g.cells.num+find(b)];
   D   = sparse(1:size(g.cellFaces,1), double(g.cellFaces(:,1)),1);
   D(:,sum(D, 1)==2) = 0;
   A   = [C, -D(:,b)]'*T(:,I1);
   
   if ~any(dF) %&& (isempty(W) || ~any(strcmpi({ W.type }, 'bhp'))),
      A(1) = A(1) + 1;
   end

   % identify internal faces
   i  = all(g.faces.neighbors ~= 0, 2);

   % Gravity contribution for each face
   cf  = g.cellFaces(:,1);
   j   = i(cf)|dF(cf);
   s   = 2*(g.faces.neighbors(cf, 1) == cellNo) - 1;
   fg  = accumarray(cf(j), grav(j).*s(j), [g.faces.num, 1]);

   
   rhs = [gg; -hh(b)];
   
   % I'm sure there's a better way to do this.
   if any(dF),
       ind        = [false(g.cells.num, 1) ; dF(b)];
       rhs        = rhs - A(:,ind)*dC;
       rhs(ind)   = dC;
       A(ind,:)   = 0;
       A(:,ind)   = 0;
       A(ind,ind) = speye(sum(ind));
   end
   
   if norm(gravity()) > 0, 
       
       
   %        rhs = rhs + T(:,I1)'*fg(cf);
   end
   
   x   = A\rhs;
   xr.cellPressure = x(1:g.cells.num);
   xr.cellFlux     = T(:, I1)*x;
   xr.faceFlux     = cellFlux2faceFlux(g, xr.cellFlux);
end 
   
%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [ff, gg, hh, grav, dF, dC] = computeRHS(g, omega, bc, src)
%Compute right-hand side contributions to pressure linear system.
%
% SYNOPSIS:
%   [f, g, h, grav, dF, dC] = computeRHS(G, omega, bc, src)
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
%   f, g, h - Direct (block) right hand side contributions as expected by
%             the linear system solvers such as 'schurComplementSymm'.
%
%   grav    - Pressure contributions from gravity.
%
%   dF, dC  - Packed Dirichlet/pressure condition structure.  The faces in
%             'dF' and the values in 'dC'.  May be used to eliminate known
%             face pressures from the linear system before calling a system
%             solver (e.g., 'schurComplementSymm').
%
% SEE ALSO:
%   addBC, addSource, pvt, computeMimeticIP, schurComplementSymm.

   ff   = zeros([size(g.cellFaces, 1), 1]);
   gg   = zeros([g.cells.num, 1]);
   hh   = zeros([g.faces.num, 1]);
   grav = grav_pressure(g, omega);

   if ~isempty(src),
      assert (max([src.cell]) <= g.cells.num);

      ss = accumarray(src.cell, src.rate);
      ii = accumarray(src.cell, 1);
      gg(ii > 0) = gg(ii > 0) + ss(ii > 0);
   end

   dF = false([g.faces.num, 1]);
   dC = [];

   if ~isempty(bc),
      assert (max([bc.face]) <= g.faces.num);
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      isDir = strcmp('pressure', bc.type);
      isNeu = strcmp('flux',     bc.type);

      [face, ix] = sort(bc.face(isDir));
      dF(face)   = true;

      % Add Dirichlet/pressure conditions to RHS.
      dC = bc.value(isDir);
      dC = dC(ix); % Order according to increasing face index.
      ff(dF(g.cellFaces(:,1))) = ff(dF(g.cellFaces(:,1))) - dC;

      % Add Neumann/flux conditions to RHS.
      % Negative sign because value is INJECTION flux.
      hh(bc.face(isNeu)) = - bc.value(isNeu);
   end

   assert (~any(dC < 0));  % Pressure conditions should always be non-neg.
end

%--------------------------------------------------------------------------

function grav = grav_pressure(g, omega)
   g_vec = gravity();   % Must be a 1-by-3 row vector for subsequent code.
   dim   = size(g.nodes.coords,2);

   assert (1 < dim && dim < 4);
   assert (all(size(g_vec) == [1,3]));

   if norm(g_vec(1 : dim)) > 0,
      cellno = rldecode(1 : g.cells.num, double(g.cells.numFaces), 2) .';
      cvec   = g.faces.centroids(g.cellFaces(:,1), :) - ...
               g.cells.centroids(cellno          , :);
      grav   = omega(cellno) .* (cvec * g_vec(1:dim)');
   else
      grav   = zeros([size(g.cellFaces,1), 1]);
   end
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['incompTPFA:', s];
end



