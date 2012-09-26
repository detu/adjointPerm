function [F, J] = blackoilJacobian(G, resSol, wellSol, pv, B, mu, tab, ...
                                   wells, src, bc)
%Residual and Jacobian of system of nonlinear equations associated with single-point upwind scheme.
%
% SYNOPSIS:
%   [F, J] = blackoilJacobian(G, resSol, wellSol, pv, B, mu, tab, ...
%                             wells, src, bc)
%
% PARAMETERS:
%   G       - Grid data structure discretizing geometric reservoir model.
%
%   resSol  - Reservoir solution structure providing fluxes, pressures and
%             fluid masses (saturations).
%
%   wellSol - Well solution structure providing fluxes and pressures for
%             all wells in the model.
%
%   pv      - Vector, size [G.cells.num, 1], of reservoir model pore
%             volumes.
%
%   B, mu   - Phase volume formation factors and fluid phase viscosities
%             for all cells.  The columns of these arrays are interpreted
%             as:
%                  1 <-> Aqua,  2 <-> Liquid,  3 <-> Vapour.
%
%   tab     - Collection of Black Oil fluid PVT and relative permeability
%             tables as defined by function 'readpvt'.
%
%   wells   - Well data structure as defined by function 'addWell'.
%
%   src, bc - Source and boundary condition structures as defined by
%             functions 'addSource' and 'addBC', respectively.
%
% RETURNS:
%   F - Anonymous function @(z,z0,dt) evaluating black oil residual vector.
%   J - Anonymous function @(z,   dt) evaluating black oil Jacobian matrix.
%
% Functions F and J are suitable for passing to function 'newtonRaphson'.
%
% SEE ALSO:
%   newtonRaphson, implicitBlackOil.

%{
#COPYRIGHT#
%}

% $Id: blackoilJacobian.m 1953 2009-03-31 10:54:12Z bska $

   assert (all(size(B) == size(mu)));

   assert (all(B (:) > 0) && ~any(isnan(B (:))));
   assert (all(mu(:) > 0) && ~any(isnan(mu(:))));

   % Compute flux matrix expansion where each flux v from cell i to cell j
   % is replaced by
   %
   %   v*spdiag([1/B(1,1), 1/B(1,2), 1/B(1,3), 1/B(2,1), ...],0)
   %
   % where the B-factors are approximated at each interface.

   [q, G] = inflow_bo(G, resSol, wellSol, tab, wells, src, bc);
   q = reshape(bsxfun(@rdivide, q, pv)', [], 1);

   % Divide by pore volume
   N = size(G,1);
   G = G / spdiags(kron(pv, ones(size(B,2),1)), 0, N,N);

   % Subtract outflow from diagonal
   G = G + spdiags(-sum(G, 1).' - min(q, 0), 0, size(G,1), size(G,2));

   % Keep positive sources in q
   q = full(max(q,0));

   % ----------------------------------------------------------------------
   % Construct return arguments Residual and Jacobian.
   % ----------------------------------------------------------------------

   function f = Residual(z, z0, dt)
      z  = reshape(z,  size(B'))';
      z0 = reshape(z0, size(B'))';
      if any(size(z) ~= size(z0)) || numel(z0) ~= numel(B),
         error('z and z0 must be same size as B');
      end
      if dt < 0,
         error('dt must be nonnegative');
      end
      f = reshape((z - z0)', [], 1) + ...
          dt * (G * reshape(flux(sat(z, B), mu)', [], 1) - q);
   end

   function jac = Jacobian(z, dt)
      z = reshape(z, size(B'))';
      if any(size(z) ~= size(B)),
         error('z and z0 must be same size as B');
      end
      if dt < 0,
         error('dt must be nonnegative');
      end
      jac = speye(numel(z)) + ...
            dt * G * dfds(sat(z, B), mu) * dsdz(z, B);
   end

   F = @Residual;
   J = @Jacobian;
end

% =========================================================================
%
% Computation of source terms and volume factors on faces. (adapted from
% inflow_bo)
%
% =========================================================================

function [q, Gm] = inflow_bo(G, resSol, wellSol, tab, wells, src, bc)

   %-----------------------------------------------------------------------
   %% Determine face pressures and masses for all faces (and wells) -------
   %
   z = resSol.z;
   p = resSol.cellPressure;

   % Wells
   nperf = 0;
   if ~isempty(wells),
      nperf = cellfun(@numel, { wells.cells }) .';
      wp    = vertcat(wellSol.pressure);

      p = [p; rldecode(wp, nperf)]; clear wp
      z = [z; resSol.z(vertcat(wells.cells),:)];  % Could be changed to well z
   end
   wi = cumsum([0; nperf]);

   % Explicit sources
   if ~isempty(src),
      c_off = G.faces.num + wi(end);
      p  = [p; resSol.cellPressure(src.cell)];
      z  = [z; resSol.z(src.cell,:)];         % Could be changed to source z
   end

   %-----------------------------------------------------------------------
   %% Phase flux matrices (pre-SPARSE format) -----------------------------
   %
   B = build_B_factors(G, tab, z, p);
   if ~isempty(bc), B = fix_B_for_bc(B, bc, G, tab, resSol); end

   % cell faces with inflow
   is_inflw = resSol.cellFlux < 0;
   fix_in   = G.cellFaces    (is_inflw, 1);
   flx_in   = resSol.cellFlux(is_inflw);

   B_inflw  = B                       (fix_in,:) ;
   N        = double(G.faces.neighbors(fix_in,:));

   %
   cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
   sgn      = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1)==cellNo)-1;
   ix       = sgn(is_inflw) > 0;
   %ix       = G.cellFaces(is_inflw,2) > 0;
   N(ix,[1,2]) = N(ix,[2,1]);

   is_int = all(N > 0, 2);

   % matrix entries for internal faces
   gi = N(is_int,2); gj = N(is_int,1);
   gs = bsxfun(@rdivide, flx_in(is_int), B_inflw(is_int,:));

   %-----------------------------------------------------------------------
   %% External sources/sinks (e.g. wells and BC's) ------------------------
   %
   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source

   % We need to treat injection into reservoir differently.  Specifically,
   % we need to take into account the fluid composition of the injected
   % fluid lest the subsequent transport solve (using e.g. 'blackoilUpwFE')
   % be thoroughly confused.
   %
   % For wells, this composition is available in the field '.comp_i' (one
   % saturation, [Aqua, Liquid, Vapour], for each *well*), while explicit
   % sources provide one three-component saturation value for each
   % injection *cell* (field 'src.sat').
   %
   if ~isempty(wells),
      % Wells as defined by 'addWell'.

      qsw = bsxfun(@rdivide, vertcat(wellSol.flux), ...
                   B(G.faces.num + (1 : wi(end)), :));

      % Injection wells have positive flux (into reservoir).
      %
      inj_w = find(cellfun(@sum, { wellSol.flux }) > 0);
      if ~isempty(inj_w),
         ix        = mcolon(wi(inj_w) + 1, wi(inj_w + 1));
         qsw(ix,:) = qsw(ix,:) .* rldecode(vertcat(wells(inj_w).compi), ...
                                           nperf(inj_w));
      end

      qi = [qi; vertcat(wells.cells)];
      qs = [qs; qsw             ];
   end

   if ~isempty(src),
      % Explicit sources defined by 'addSource'.

      qss = bsxfun(@rdivide, src.rate, B(c_off + 1 : end, :));

      % Injection sources have positive rate into reservoir.
      %
      in  = find(src.rate > 0);
      if ~isempty(in),
         qss(in,:) = qss(in,:) .* src.sat(in,:);
      end

      qi = [qi; src.cell];
      qs = [qs; qss     ];
   end

   if ~isempty(bc),
      % Contributions from boundary conditions as defined by 'addBC'.

      qsbc  = zeros([G.faces.num, 1]);
      dirF  = false([G.faces.num, 1]);

      isDir = strcmp('pressure', bc.type);
      isNeu = strcmp('flux',     bc.type);

      dirF(bc.face(isDir)) = true;

      % Depends on Dirichlet faces being outer faces.
      qsbc(dirF)           = -resSol.cellFlux(dirF(G.cellFaces(:,1)));
      %qsbc(dirF)           = -resSol.cellFlux(logical(S.D * dirF));
      qsbc(bc.face(isNeu)) = bc.value(isNeu);

      % Injection BC's have positive rate into reservoir.
      %
      is_inj   = find(qsbc > 0);
      inj_comp = bc.sat(ismember(bc.face, is_inj), :);

      qsbc = bsxfun(@rdivide, qsbc, B(1:G.faces.num, :));

      if ~isempty(is_inj),
         qsbc(is_inj,:) = qsbc(is_inj,:) .* inj_comp;
      end

      is_outer = any(G.faces.neighbors == 0, 2);

      qi = [qi; sum(G.faces.neighbors(is_outer,:), 2)];
      qs = [qs; qsbc(is_outer,:)];
   end

   %-----------------------------------------------------------------------
   %% Assemble final phase flux and source contributions in SPARSE format -
   %
   nc   = G.cells.num;     % Number of cells
   np   = size(B,2);       % Number of phases
   nsrc = numel(qi);       % Number of cells containing explicit sources

   % The matrix 'Gm' collects the mass flux for all phases in each cell
   % into blocks of contiguous rows and columns.  Specifically, the mass
   % fluxes (of all the 'np' phase) from cell 'j' and into cell 'i'
   % constitute the entries
   %
   %    Gm((i-1)*np + (1 : np), (j-1)*np + (1 : np))
   %
   % of the mass flux matrix.
   %
   % Assume 'x' is a column vector of (positive) integer values in the
   % range [1 : G.cells.num].  Then
   %
   %    I = BSXFUN(@PLUS, (1 : np).', np * (x - 1).')
   %
   % is a matrix of size [np, NUMEL(x)] the column 'c' of which contains
   % the linear indices
   %
   %    (x(c) - 1)*np + (1 : np)
   %
   % Using the 'ix' indexing function, the single call to SPARSE
   % assembles the mass flux matrix.
   %
   ix = @(x) reshape(bsxfun(@plus, (1 : np).', np * (x - 1).'), [], 1);
   Gm = sparse (ix(gi), ix(gj), reshape(gs', [], 1), nc*np, nc*np);

   qj = reshape(repmat(1:np, [nsrc, 1]), [], 1);
   if ~isempty(qi),
      q = sparse(reshape(qi(:,ones([1,np])),[],1), qj, qs(:), nc, np);
   else
      q = sparse(nc,np);
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


function B = build_B_factors(G, tab, z, p)
   [B, B] = evalpvt(tab, z, p ./ 1e5);
   B = [zeros([1, size(z,2)]); B];

   B1 = B(G.faces.neighbors(:,1) + 1, :);
   B2 = B(G.faces.neighbors(:,2) + 1, :);
   Bf = bsxfun(@rdivide, ...
               sum(cat(3, B1, B2), 3), ...
               sum(G.faces.neighbors > 0, 2));
   B  = [Bf; B(G.cells.num + 2 : end, :)];
end

%--------------------------------------------------------------------------

function B = fix_B_for_bc(B, bc, G, tab, resSol)
   faces = [bc.face];
   N = G.faces.neighbors(faces,:);  assert (all(xor(N(:,1)==0, N(:,2)==0)));

   fp = resSol.facePressure(faces);
   z  = resSol.z(sum(N,2), :);
   [cf, Bf] = evalpvt(tab, z, fp ./ 1e5);
   B(faces,:) = Bf;
end


% =========================================================================
%
% Computation of functions and derivatives needed by Residual and Jacobian
%
% =========================================================================


% -------------------------------------------------------------------------
function f = flux(s, mu)
kr      = getRelPerm(s);
L       = bsxfun(@rdivide, kr, mu);
f       = bsxfun(@rdivide, L, sum(L,2));
end

% -------------------------------------------------------------------------

function df = dfds(s, mu)
% df = dfds(s, fluid) - compute jacobian of fractional flow function
%
% Since sum(f, 2)==1 and sum(s, 2)==1 by definition, the rank of df is
% size(s, 2)-1, i.e., on less than th number of phases present.
%
%             kr_i/mu_i          L_i
%  f_i  = ----------------- = -----------
%          \sum_j kr_j/mu_j      L_T
%
%  df_i   dL_i/ds_j      L_i
%  ---- = --------- -  -------- * dL_T/ds_j   (*)
%  ds_j     L_T          L_T^2

% Compute relative permeability and derivatives
[kr, dkr] = getRelPerm(s);

% Mobility
L         = bsxfun(@rdivide, kr, mu);

% Inverse of cell total mobility diagonal matrix
iLt       = spdiags(1 ./ kron(sum(L, 2), ones(size(s,2), 1)), 0, ...
                    numel(mu), numel(mu));

% Inverse of phase viscosity diagonal matrix
assert (size(mu, 1) == size(s, 1));
imu       = spdiags(1./reshape(mu', [], 1), 0, numel(mu), numel(mu));


df1       = imu*dkr;

% derivative of total mobility
dLt       = sum(df1, 2);

% First term of (*)
df1       = iLt*df1;

% Second term of (*)
df2       = -blockOuterProduct(size(s,2), iLt*iLt*reshape(L', [], 1), dLt);

df = df1 + df2;
end

% -------------------------------------------------------------------------
function s  = sat(z, B)
% s = sat(z, fluid) - compute saturation, given mass distribution z
%
%           z_i*B_i
%  s_i = ------------
%        \sum z_iB_i
%

u = bsxfun(@times, z, B);
s = bsxfun(@rdivide, u, sum(u,2));
end

% -------------------------------------------------------------------------
function ds = dsdz(z, B)
% ds = dsdz(z, fluid) - compute jacobian of saturation wrt. mass
% distribution z.
%
% Since sum(s, 2)==1 by definition, the rank of dsdz is size(s,2)-1.
%
%         z_i*B_i
% s_i = ------------------
%        \sum_j z_j*B_j
%
%
% ds_i        B_i               z_i*B_i  B_j
% ---- = ------------------- - ---------------------
% dz_j    \sum_j z_j*B_j         (\sum_j z_j*B_j)^2

%

% Phase volumes
u   = bsxfun(@times, z, B);


ds1 =  spdiags(kron(1./sum(u,2), ones(size(B,2),1)).*reshape(B',[], 1), ...
               0, numel(u), numel(u));
ds2 = -blockOuterProduct(size(z, 2), u', bsxfun(@rdivide, B, sum(u,2).^2)');
ds  = ds1 + ds2;
end

% -------------------------------------------------------------------------
function A = blockOuterProduct(blocksize, u, v)
nblocks = numel(u)/blocksize;
i       = kron(1:nblocks, ones([blocksize,1]));
j       = 1:nblocks*blocksize;
U       = sparse(i, j, u);
V       = sparse(i, j, v);
A       = U'*V;
end
