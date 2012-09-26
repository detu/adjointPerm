function [F, J, rhov] = blackoilJacobian2(G, resSol, wellSol, pv, fluid, ...
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

% $Date: 2009-10-05 18:06:36 +0200 (ma, 05 okt 2009) $
% $Revision: 2951 $

   
   % Compute flux matrix expansion where each flux v from cell i to cell j
   % is replaced by
   %
   %   v*spdiag([1/B(1,1), 1/B(1,2), 1/B(1,3), 1/B(2,1), ...],0)
   %
   % where the B-factors are approximated at each interface.

   nc = size(resSol.z, 1);
   np = size(resSol.z, 2);
   
   [q, rhov] = inflow_bo_new(G, resSol, wellSol, fluid, wells, src, bc);
   q = reshape(bsxfun(@rdivide, q, pv)', [], 1);

   % Divide by pore volume
   N = size(rhov,1);
   rhov = spdiags(kron(pv, ones([np,1])), 0, N, N) \ rhov;

   % Keep positive sources in q
   q = full(q);

   % cell pressure
   p = resSol.cellPressure;
   
   % ----------------------------------------------------------------------
   % Construct return arguments Residual and Jacobian.
   % ----------------------------------------------------------------------
   % ----------------------------------------------------------------------
   function f = flux(s, mu)
      f = bsxfun(@rdivide, fluid.relperm(s), mu);
      f = bsxfun(@rdivide, f, f * ones([size(f,2), 1]));
   end

   function f = Residual(z, z0, dt)
      % Assume z (np*nc-by-1) is ordered per phase, then per cell.
      z  = reshape(z,  [np, nc])';
      z0 = reshape(z0, [np, nc])';
      if any(size(z) ~= size(z0)) || numel(z0) ~= np*nc,
         error('z and z0 must be same size');
      end
      if dt < 0,
         error('dt must be nonnegative');
      end
      [c, rho, mu, u, R] = fluid.pvt(p, z); %#ok
      
      f = reshape((z - z0)', [], 1) + ...
          dt * (rhov * reshape(flux(sat(u), mu)', [], 1) - q);
   end

   function jac = Jacobian(z, dt)
      z = reshape(z, [np, nc])';
      if any(size(z) ~= [nc, np]),
         error('z and z0 must be same size');
      end
      if dt < 0,
         error('dt must be nonnegative');
      end
      [c, rho, mu, u, R] = fluid.pvt(p, z);
      rhoS = fluid.surfaceDensity;
   
      jac = speye(numel(z)) + ...
            dt * rhov * dfds(sat(u), mu, @fluid.relperm) * dsdu(u) *...
            dudz(rho, rhoS, R);
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
function B = computeFaceDensity(G, fluid, resSol, wellSol, wells, src)
   % Use upstream definition of 'z' for internal faces.
   i  = all(G.faces.neighbors > 0, 2);
   N  = G.faces.neighbors(i,:);
   fn = resSol.faceFlux(i) < 0;

   N(fn, :)  = N(fn, [2,1]);
   up_ix     = zeros([G.faces.num, 1]);
   up_ix( i) = N(:, 1);
   up_ix(~i) = sum(G.faces.neighbors(~i,:), 2);

   p = resSol.facePressure;
   z = resSol.z(up_ix, :);

   % Wells
   if ~isempty(wells),
      nperf = cellfun(@numel, { wells.cells }) .';
      wp    = vertcat(wellSol.pressure);

      p = [p;rldecode(wp, nperf)]; clear wp
      z = [z;resSol.z(vertcat(wells.cells),:)];  % Could be changed to well z
   end

   % Explicit sources
   if ~isempty(src),
      p  = [p; resSol.cellPressure(src.cell)];
      z  = [z; resSol.z(src.cell,:)];         % Could be changed to source z
   end

   [rhof, rhof] = fluid.pvt(p, z);
   B = bsxfun(@rdivide, fluid.surfaceDensity, rhof);
end

%--------------------------------------------------------------------------


function [q, rhov] = inflow_bo_new(G, resSol, wellSol, fluid, wells, src, bc)

   %-----------------------------------------------------------------------
   %% Determine face pressures and masses for all faces (and wells) -------
   B = computeFaceDensity(G, fluid, resSol, wellSol, wells, src);
   
   flux     = bsxfun(@rdivide, resSol.faceFlux, B(1:G.faces.num, :));
   N        = double(G.faces.neighbors);
   %cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
   %sgn      = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1)==cellNo)-1;
   
   gi = 1+N(:);
   gj = 1+reshape(N(:,[2,1]), [], 1);
   gs = [flux; -flux];
   
   % Find all inflow fluxes on boundary
   qi = [];
   qs = [];
   if ~isempty(bc),
      sgn = 1-2*(N(bc.face, 1)==0);
      bcflux = -bsxfun(@times, flux(bc.face,:), sgn);
      inflow = resSol.faceFlux(bc.face).*sgn < 0;
      qi = sum(N(bc.face(inflow), :), 2);
      qs = bcflux(inflow, :);
      
      if any(inflow),
         qs = qs.*bc.sat(inflow, :);
      end
   end
   
   % Find outflow. Mover outflow to diagonal
   outflow     = gs(:,1)>0;
   gj(outflow) = gi(outflow);
   
   % Keep internal inflow fluxes and all outflow fluxes
   ind = gi~=1 & gj(:)~=1;
   gi = gi(ind)-1; 
   gj = gj(ind)-1;
   gs = gs(ind,:);
   
   
   %% External sources/sinks (e.g. wells and BC's) ------------------------
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
   wi = 0;
   if ~isempty(wells),
      % Wells as defined by 'addWell'.
      nperf = cellfun(@numel, { wells.cells }) .';
      wi = cumsum([0;  nperf]);
      qsw = bsxfun(@rdivide, vertcat(wellSol.flux), ...
                   B(G.faces.num + (1 : wi(end)), :));
      %qsw = repmat(vertcat(wellSol.flux), [1,size(B, 2)]);
      % Injection wells have positive flux (into reservoir).
      %
      inj_w = find(cellfun(@sum, { wellSol.flux }) > 0);
      zperf = rldecode(vertcat(wells.compi), nperf);
      pperf = rldecode(vertcat(wellSol.pressure), nperf);
      vperf = repmat(vertcat(wellSol.flux), [1,3]);
      if ~isempty(inj_w),
         ix        = mcolon(wi(inj_w) + 1, wi(inj_w + 1));
         
         [u, u, u, u] = fluid.pvt(pperf(ix), zperf(ix,:));
         alpha        = sum(u, 2);
         qsw(ix,:) = vperf(ix,:) .* bsxfun(@rdivide, zperf(ix, :), alpha);
      end

      qi = [qi; vertcat(wells.cells)];
      qs = [qs; qsw             ];
   end

   if ~isempty(src),
      % Explicit sources defined by 'addSource'.

      qss = bsxfun(@rdivide, src.rate, B(G.faces.num + wi(end) + 1 : end, :));

      % Injection sources have positive rate into reservoir.
      %
      in  = find(src.rate > 0);
      if ~isempty(in),
         qss(in,:) = qss(in,:) .* src.sat(in,:);
      end

      qi = [qi; src.cell];
      qs = [qs; qss     ];
   end
   
   % hack 
   ind = qs(:,1) < 0;
   if any(qs(ind,:)>0, 2), 
      error ('For some reason, the phase source terms differ in sign.');
   end
   gs  = [gs; -qs(ind, :)];
   gi  = [gi; qi(ind)];
   gj  = [gj; qi(ind)];
   qi(ind)   = [];
   qs(ind,:) = [];
   
   %-----------------------------------------------------------------------
   %% Assemble final phase flux and source contributions in SPARSE format -
   %
   nc   = G.cells.num;     % Number of cells
   np   = size(B,2);       % Number of phases
   nsrc = numel(qi);       % Number of cells containing explicit sources

   % The matrix 'rhov' collects the mass flux for all phases in each cell
   % into blocks of contiguous rows and columns.  Specifically, the mass
   % fluxes (of all the 'np' phase) from cell 'j' and into cell 'i'
   % constitute the entries
   %
   %    rhov((i-1)*np + (1 : np), (j-1)*np + (1 : np))
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
   rhov = sparse (ix(gi), ix(gj), reshape(gs', [], 1), nc*np, nc*np);

   qj = reshape(repmat(1:np, [nsrc, 1]), [], 1);
   if ~isempty(qi),
      q = sparse(reshape(qi(:,ones([1,np])),[],1), qj, qs(:), nc, np);
   else
      q = sparse(nc,np);
   end
end

% =========================================================================
%
% Computation of functions and derivatives needed by Residual and Jacobian
%
% =========================================================================

% -------------------------------------------------------------------------

function df = dfds(s, mu, relperm)
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
[kr, dkr] = relperm(s);

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
function s  = sat(u)
% s = sat(z, fluid) - compute saturation, given mass distribution z
%
%            u_i
%  s_i = ------------
%          \sum u_i

   s = bsxfun(@rdivide, u, u * ones([size(u,2), 1]));
end

% -------------------------------------------------------------------------
function ds = dsdu(u)
% ds = dsdz(z, fluid) - compute jacobian of saturation wrt. mass
% distribution z.
%
% Since sum(s, 2)==1 by definition, the rank of dsdz is size(s,2)-1.
%
%         u_i
% s_i = ------------------
%        \sum_j u_j
%
%
% ds_i        del_ij          u_i  
% ---- = ------------- - --------------
% du_j    \sum_k u_k     (\sum_k u_k)^2

   nc = size(u,1);
   np = size(u,2);
   
   ds1 =  spdiags(kron(1./sum(u,2), ones(np,1)), 0, np*nc, np*nc);
   ds2 = -blockOuterProduct(np, u', repmat(1./sum(u,2).^2, [1,np])' );
   ds  = ds1 + ds2;
end

% -------------------------------------------------------------------------
function du = dudz(rho, rhoS, R)
   % du = dudz(B, R) - compute change in phase volume when mass is altered
   %
   % u = inv(B)*R*z.  Let dudz be approximated by inv(B)*R.
   
   
   % nei:  u = B*inv(R)*z, dudz = B*inv(R) + (d/dz(B)*inv(R) + B*d/dz(inv(R)))*diag(z)
   % we ignore the second term for now.
   
   nc = size(rho,1);
   np = size(rho,2);
   
   B  = bsxfun(@rdivide, rhoS, rho);
   d  = 1-R(:,1).*R(:,2);
   i  = bsxfun(@plus, np*(0:nc-1), [1,2,3,1,2,2,1]');
   j  = bsxfun(@plus, np*(0:nc-1), [1,2,3,2,1,3,3]');
   v  = [B, -B(:,[2,3,3]).*R(:,[3,2,1]), B(:,2).*R(:,1).*R(:,3)];
   v  = bsxfun(@rdivide, v, d) .';
   du = sparse(i, j, v, nc*np, nc*np);
   
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


%{
function [q, rhov] = inflow_bo(G, resSol, wellSol, fluid, wells, src, bc)

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
   B = build_B_factors(G, fluid, z, p);
   if ~isempty(bc), B = fix_B_for_bc(B, bc, G, fluid, resSol); end

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
   
   N(ix,[1,2]) = N(ix,[2,1]);

   is_int = all(N > 0, 2);

   % matrix entries for internal faces
   gi = N(is_int,2); gj = N(is_int,1);
   gs = bsxfun(@rdivide, flx_in(is_int), B_inflw(is_int,:));

   
   % nytt forsøk
   massflux = bsxfun(@rdivide, resSol.faceFlux, B);
   
   
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
      %qsw = repmat(vertcat(wellSol.flux), [1,size(B, 2)]);
      % Injection wells have positive flux (into reservoir).
      %
      inj_w = find(cellfun(@sum, { wellSol.flux }) > 0);
      zperf = rldecode(vertcat(wells.compi), nperf);
      pperf = rldecode(vertcat(wellSol.pressure), nperf);
      vperf = repmat(vertcat(wellSol.flux), [1,3]);
      if ~isempty(inj_w),
         ix        = mcolon(wi(inj_w) + 1, wi(inj_w + 1));
         
         [u, u, u, u] = fluid.pvt(pperf(ix), zperf(ix,:));
         alpha        = sum(u, 2);
         qsw(ix,:) = vperf(ix,:) .* bsxfun(@rdivide, zperf(ix, :), alpha);
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

   % The matrix 'rhov' collects the mass flux for all phases in each cell
   % into blocks of contiguous rows and columns.  Specifically, the mass
   % fluxes (of all the 'np' phase) from cell 'j' and into cell 'i'
   % constitute the entries
   %
   %    rhov((i-1)*np + (1 : np), (j-1)*np + (1 : np))
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
   rhov = sparse (ix(gi), ix(gj), reshape(gs', [], 1), nc*np, nc*np);

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


function B = build_B_factors(G, fluid, z, p)
   [rho, rho] = fluid.pvt(p, z);
   rhoS = fluid.surfaceDensity;

   rho  = [zeros([1, size(z,2)]); rho];
   rho1 = rho(G.faces.neighbors(:,1) + 1, :);
   rho2 = rho(G.faces.neighbors(:,2) + 1, :);
   rhof = bsxfun(@rdivide, ...
                 sum(cat(3, rho1, rho2), 3), ...
                 sum(G.faces.neighbors > 0, 2));
   Bf = bsxfun(@rdivide, rhoS, rhof);
   B  = [Bf; bsxfun(@rdivide, rhoS, rho(G.cells.num + 2 : end, :))];
end

%--------------------------------------------------------------------------

function B = fix_B_for_bc(B, bc, G, fluid, resSol)
   faces = [bc.face];
   N = G.faces.neighbors(faces,:);  assert (all(xor(N(:,1)==0, N(:,2)==0)));

   fp = resSol.facePressure(faces);
   z  = resSol.z(sum(N,2), :);
   
   [rho, rho] = fluid.pvt(fp, z);
   Bf         = bsxfun(@rdivide, fluid.surfaceDensity, rho);
   
   B(faces,:) = Bf;
end
%}
