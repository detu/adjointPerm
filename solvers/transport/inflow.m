function [gm, sources, porvol] = inflow(G, rock, S, W, resSol, wellSol, ...
                                        varargin)
%Determine upwind fluxes (and sources) from flow calculation results.
%
% SYNOPSIS:
%   [gm, sources, porvol] = inflow(G, rock, S, W, resSol, wellSol)
%   [gm, sources, porvol] = inflow(G, rock, S, W, resSol, wellSol,
%                                  'pn1', pv1, ...)
% DESCRIPTION:
%
% PARAMETERS:
%   G       - Grid data structure as defined by 'grid_structure'.
%
%   rock    - Rock data structure.  Must contain the field 'rock.poro'.
%
%   S       - Aggregate system structure, particularly containing
%             (non-homogeneous) Neumann (flux) boundary conditions in field
%             S.RHS.h .
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
% RETURNS:
%   gm      - In-flow connection matrix suitable for passing to function
%             twophaseUpwBE.
%
%   sources - Aggregate source term suitable for passing to function
%             twophaseUpwBE.
%
%   porvol  - Vector of size G.cells.num-by-1 of cell pore volumes.
%
% TODO:
%           - Handling of sources is cut and adapted from inflow_bo,
%             would be better to make a utility function.
%
% SEE ALSO:
%   solveBlackOilWellSystem, solveWellSystemMS, twophaseUpwBE.

%{
#COPYRIGHT#
%}

% $Id: inflow.m 1953 2009-03-31 10:54:12Z bska $

assert(~isempty(rock));

opt = struct('src', [], 'bc', []);
opt = merge_options(opt, varargin{:});
src = opt.src;
bc  = opt.bc;

% resSol.cellFlux is amount of outflow for all half-faces on all cells.
% Negative means inflow.
%
is_inflow    = resSol.cellFlux < 0;

inflow_faces = G.cellFaces(is_inflow, 1);
inflow_flux  = resSol.cellFlux(is_inflow);

N            = double(G.faces.neighbors(inflow_faces, :));

cellNo       = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
sgn          = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1)==cellNo)-1;

%ind         = G.cellFaces(is_inflow, 2) > 0;
ind          = sgn(is_inflow) > 0;

N(ind,[1,2]) = N(ind,[2,1]);

% In distributing the flow pattern, only consider internal faces when
% constructing the connection matrix.  Flux boundary conditions, if any,
% are considered sources and handled below.
%
int_face = prod(N, 2) > 0;
int_conn = N(int_face, :);
int_flux = inflow_flux(int_face);

nc = G.cells.num;
gm = sparse(int_conn(:,2), int_conn(:,1), int_flux, nc, nc);

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
