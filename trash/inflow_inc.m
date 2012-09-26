function [gravity, flux, sources, porvol, dt] = ...
    inflow_inc(G, rock, fluid, S, W, resSol, wellSol, varargin)
% INFLOW_INC -- Compute input to transport solver from solution to pressure
% equation. If specified, also determine timestep for explicit solver.
%
% matrix.
% SYNOPSIS:
%   [gravity, flux, sources, porvol, dt] =
%        inflow_inc(G, rock, S, W, resSol, wellSol)
%   [gravity, flux, sources, porvol, dt] =
%        inflow_inc(G, rock, S, W, resSol, wellSol, 'pn1', pv1, ...)
% DESCRIPTION:
%
% PARAMETERS:
%   G       - grid_structure grid structure.
%
%   rock    - Rock data structure.  Must contain the field 'rock.poro'.
%
%   S       - Aggregate system structure.
%
%   W       - Well data structure as defined by function addWell.
%
%   resSol  - Reservoir solution structure.
%
%   wellSol - Well solution structure.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
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
%               - OnlyGrav  --  Whether or not to only compute gravity
%                               contribution (for use in gravity
%                               splitting). Default value: false.
%
%               - ComputeDt --  Whether or not to compute timestep.
%                               Default value: false.
%
%               - max_dt    --  Maximal timestep. Default value = 100;
%
%               - dt_factor --  Factor to multiply dt by.
%                               Default value = 0.5.
%
% RETURNS:
%   flux    - In-flow connection matrix suitable for passing to function
%             twophaseUpwFEGrav. Size G.cells.num-by-(G.faces.num-boundary
%             faces)
%
%   gravity - Matrix with gravity contribution for each face. Same size as
%             flux.
%
%   sources - Aggregate source term suitable for passing to function
%             twophaseUpwFEGrav.
%
%   porvol  - Vector of size G.cells.num-by-1 of cell pore volumes.
%
%   dt      - Timestep for use with twophaseUpwFEGrav
%
% SEE ALSO:
%   solveMimeticWellSystem, twophaseUpwFEGrav, twophaseUpwBE,
%   twophaseUpwBEGrav.

% TODO:
%        - use mimetic permeability instead of harmonic perm
%        - add support for bc (dir of flow might change when g~=0)
%        - Handling of sources is cut and adapted from inflow_bo,
%          would be better to make a utility function.

% $Id: inflow_inc.m 1829 2009-03-23 14:01:07Z bska $

opt = struct('OnlyGrav', false, ...
             'ComputeDt', false, ...
             'max_dt' , 100, ...
             'dt_factor', 0.5, ...
             'src', [], ...
             'bc', []);

opt       = merge_options(opt, varargin{:});
computeDt = opt.ComputeDt;
max_dt    = opt.max_dt;
dt_factor = opt.dt_factor;
onlyGrav  = opt.OnlyGrav;
src       = opt.src;
bc        = opt.bc;

cellNo = rldecode((1:G.cells.num)', double(G.cells.numFaces));
%to remove boundary faces - temporary until bc are included:
intFInx = prod(double(G.faces.neighbors),2)~=0;
cIntFInx = logical(S.D * intFInx);

[i, trash] = find(S.D(:,intFInx)'); %#ok
%-----------------------------------------------------------------------
%% Assemble flux matrix ------------------------------------------------
%
if onlyGrav,
   flux = sparse(G.cells.num, nnz(intFInx));
else
   flux = sparse(cellNo(cIntFInx), i, resSol.cellFlux(cIntFInx), ...
                 G.cells.num, nnz(intFInx));
end

%-----------------------------------------------------------------------
%% Assemble gravity matrix ---------------------------------------------
%
g = S.constants.gravity*(fluid.rho(1)-fluid.rho(2))*S.constants.FCON;

if g ~=0
   %spread gravity to cellfaces: cellFInx * g_const .* cellFace_normal
   Grav = (S.D*g*G.faces.normals(:,3)).*double(G.cellFaces(:,2));
   gravity = sparse(cellNo(cIntFInx), i, Grav(cIntFInx), ...
                    G.cells.num, nnz(intFInx));

   %use harmonic permeability for edges
   perm = (rock.perm(G.faces.neighbors(intFInx,1),end).* ...
           rock.perm(G.faces.neighbors(intFInx,2),end))./...
          (rock.perm(G.faces.neighbors(intFInx,1),end)+...
           rock.perm(G.faces.neighbors(intFInx,2),end));
   gravity = bsxfun(@times, perm', gravity);
else
   gravity = sparse(G.cells.num, nnz(intFInx));
end
%-----------------------------------------------------------------------
%% External sources/sinks (e.g. wells and BC's) ------------------------
%
qi = [];  % Cells to which sources are connected
qs = [];  % Actual strength of source

% We need to treat injection into reservoir differently.  Specifically,
% we need to take into account the fluid composition of the injected
% fluid lest the subsequent transport solve (using e.g. 'twophaseUpwFE')
% be thoroughly confused.
%
% For wells, this composition is available in the field '.comp_i' (one
% saturation, [Aqua, Liquid, Vapor], for each *well*), while explicit
% sources provide one three-component saturation value for each
% injection *cell* (field 'src.sat').
%
if ~isempty(W),
   % Wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';
   wi = cumsum([0; nperf]);

   qsw = vertcat(wellSol.flux);
   % Injection wells have positive flux (into reservoir).
   %
   inj_w = find(cellfun(@sum, { wellSol.flux }) > 0);
   if ~isempty(inj_w),
      ix      = mcolon(wi(inj_w) + 1, wi(inj_w + 1));
      comp    = rldecode(vertcat(W(inj_w).compi), nperf(inj_w));
      qsw(ix) = qsw(ix) .* comp(:,1);
   end

   qi = [qi; vertcat(W.cells)];
   qs = [qs; qsw             ];
end

if ~isempty(src),
   % Explicit sources defined by (e.g.) 'addSource'.

   qss = src.rate;
   % Injection sources have positive rate into reservoir.
   %
   in  = find(src.rate > 0);
   if ~isempty(in),
      qss(in) = qss(in) .* src.sat(in,1);
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
   inj_comp = bc.sat(ismember(bc.face, is_inj), 1);

   if ~isempty(is_inj),
      qsbc(is_inj) = qsbc(is_inj) .* inj_comp;
   end

   is_outer = any(G.faces.neighbors == 0, 2);

   qi = [qi; sum(G.faces.neighbors(is_outer), 2)];
   qs = [qs; qsbc(is_outer)];
end

%-----------------------------------------------------------------------
%% Assemble final phase flux and source contributions in SPARSE format -
%
sources = sparse(qi, 1, qs, G.cells.num, 1);
porvol  = poreVolume(G, rock);


%-----------------------------------------------------------------------
%%  COMPUTE CFL-number/dt ----------------------------------------------
%
dt = max_dt;
if computeDt
    % NB: only consider gravitation
    if g ~= 0
        tmp = linspace(0,1,100);
        f  =  @(fluid,resSol) (prod(fluid.kr(resSol),2)./prod(fluid.mu,2))...
            ./ fluid.Lt(resSol);
        grav_in = gravity;
        tmpSol.s(:,1) = tmp;
        grav_in(gravity>0) = 0;
        grav_in = abs(sum(grav_in,2));
        is_grav = grav_in ~=0;
        dt = dt_factor*min(porvol(is_grav)./grav_in(is_grav)) ./ ...
                       max(abs(diff(f(fluid,tmpSol))./diff(tmp)'));
    end

    % -- compute dt dependent of saturation -- %
    s = resSol.s(:,1);
    dt = min(dt,1e99);
    %interior
    is_int = all(G.faces.neighbors > 0, 2);
    N_int  = G.faces.neighbors(is_int, :) ;
    fac = dpdS(s(N_int(:,1),:), s(N_int(:,2),:), fluid);
    dt  = min(dt, min( min(porvol(N_int),[], 2)./(fac.*  abs(resSol.faceFlux(is_int)) ) ));
    %sources
    dt = min(dt,min(porvol./(abs(sources))));
    dt = min(max_dt,dt*dt_factor);
end
    function fac = dpdS(sat1,sat2,fluid)
        f1 = fractionalFlow(sat1,fluid);
        f2 = fractionalFlow(sat2,fluid);
        fac = (f1-f2);
        den = (sat1-sat2);
        is_ok = all((den < 0) + (den > 0), 2);
        fac(is_ok,:) = fac(is_ok,:) ./ den(is_ok,:);
        fac = max(abs(fac), [], 2);
    end
    function f = fractionalFlow(sat, fluid)
        sol.s = sat;
        kr = fluid.kr(sol);
        f = (kr(:,1)./fluid.mu(1))./fluid.Lt(sol);
    end
end
