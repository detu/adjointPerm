function [gm, q, B] = inflow_bo(G, S, W, resSol, wellSol, tab, varargin)
%Determine upwind Black-Oil fluxes (and sources) from flow calucations.
%
% SYNOPSIS:
%   [gm, q] = inflow_bo(G, S, W, resSol, wellSol, tab)
%   [gm, q] = inflow_bo(G, S, W, resSol, wellSol, tab, 'pn1', pv1, ...)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G       - Grid data structure as defined by 'grid_structure'.
%
%   S       - Aggregate system structure, particularly containing
%             (non-homogeneous) Neumann (flux) boundary conditions in field
%             S.RHS.h_bc.
%
%   W       - Well data structure as defined by function addWell.
%
%   resSol  - Reservoir solution structure.
%
%   wellSol - Well solution structure.
%
%   tab     - Collection of Black Oil fluid PVT and relative permeability
%             tables as defined by function 'readpvt'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - bc  -- Boundary condtion structure as defined by function
%                        'addBC'.  This structure accounts for all external
%                        boundary contributions to the reservoir flow.
%                        Default value: bc = [] meaning all external
%                        no-flow (homogeneous Neumann) conditions.
%
%               - src -- Explicit source contributions as defined by
%                        function 'addSource'.
%                        Default value: src = [] meaning no explicit
%                        sources exist in the model.
%
% RETURNS:
%   gm      - In-flow connection matrix suitable for passing to function
%             blackoilUpwBE.
%
%   q       - Aggregate source term suitable for passing to function
%             blackoilUpwFE.
%
% SEE ALSO:
%   solveBlackOilWellSystem, solveWellSystemMS, blackoilUpwBE.

%{
#COPYRIGHT#
%}

% $Id: inflow_bo.m 1953 2009-03-31 10:54:12Z bska $

   opt = struct('src', [], 'bc', []);
   opt = merge_options(opt, varargin{:});

   src = opt.src;
   bc  = opt.bc;

   %-----------------------------------------------------------------------
   %% Determine face pressures and masses for all faces (and wells) -------
   %
   z = resSol.z;
   p = resSol.cellPressure;

   % Wells
   nperf = 0;
   if ~isempty(W),
      nperf = cellfun(@numel, { W.cells }) .';
      wp    = vertcat(wellSol.pressure);

      p = [p; rldecode(wp, nperf)]; clear wp
      z = [z; resSol.z(vertcat(W.cells),:)];
   end
   wi = cumsum([0; nperf]);

   % Explicit sources
   if ~isempty(src),
      c_off = G.faces.num + wi(end);
      p  = [p; resSol.cellPressure(src.cell)];
      z  = [z; resSol.z(src.cell,:)];
   end

   %-----------------------------------------------------------------------
   %% Phase flux matrices (pre-SPARSE format) -----------------------------
   %
   B = build_B_factors(G, tab, z, p);
   if ~isempty(bc), B = fix_B_for_bc(B, bc, G, tab, resSol); end

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
   if ~isempty(W),
      % Wells as defined by 'addWell'.

      qsw = bsxfun(@rdivide, vertcat(wellSol.flux), ...
                   B(G.faces.num + (1 : wi(end)), :));

      % Injection wells have positive flux (into reservoir).
      %
      inj_w = find(cellfun(@sum, { wellSol.flux }) > 0);
      if ~isempty(inj_w),
         ix        = mcolon(wi(inj_w) + 1, wi(inj_w + 1));
         qsw(ix,:) = qsw(ix,:) .* rldecode(vertcat(W(inj_w).compi), ...
                                           nperf(inj_w));
      end

      qi = [qi; vertcat(W.cells)];
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
      qsbc(dirF)           = -resSol.cellFlux(logical(S.D * dirF));
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

   gm = cell ([np     , 1]);
   qj = zeros([np*nsrc, 1]);

   for k = 1 : np,
      gm{k}  = sparse(gi, gj, gs(:,k), nc, nc);
      qj((k-1)*nsrc + 1 : k*nsrc) = k;
   end

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
