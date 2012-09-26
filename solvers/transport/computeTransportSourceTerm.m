function q = computeTransportSourceTerm(resSol, wellSol, G, wells, src, bc)
%Compute source terms for transport
%
% SYNOPSIS:
%   resSol = getSources(resSol, wellSol, G, wells, src, bc)
%
% DESCRIPTION:
%
% REQUIRED PARAMETERS:
%   resSol  - 
%
%   wellSol - 
%
%   G       - 
%
%   wells   - 
%
%   src     - 
%
%   bc      - 
%
% RETURNS:
%   q       - 
%
% SEE ALSO:
%   twophaseJacobian.

% TODO:
%   - implement gravity effects for pressure boundary and wells

%{
#COPYRIGHT#
%}

% $Date: 2009-10-12 15:33:42 +0200 (ma, 12 okt 2009) $
% $Revision: 2980 $

%% External sources/sinks (e.g. wells and BC's) ------------------------
   %
   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source term (in m^3/s).

   if ~isempty(wells),
      [i, s] = contrib_wells(wells, wellSol);
      qi = [qi; i];
      qs = [qs; s];
   end

   if ~isempty(src), assert (~isempty(src.sat))
      [i, s] = contrib_src(src);
      qi = [qi; i];
      qs = [qs; s];
   end
   
   if ~isempty(bc), assert (~isempty(bc.sat))
      is_int = all(double(G.faces.neighbors) > 0, 2);
      [i, s] = contrib_bc(G, resSol, bc, is_int);
      qi = [qi; i];
      qs = [qs; s];
   end

   %-----------------------------------------------------------------------
   %% Assemble final phase flux and source contributions in SPARSE format -
   %
   q  = sparse(qi, 1, qs, G.cells.num, 1);
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_wells(W, wellSol)
   % Wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';

   qi = vertcat(W.cells);
   qs = vertcat(wellSol.flux);

   % Injection perforations have positive flux (into reservoir).
   %
   comp      = rldecode(vertcat(W.compi), nperf);
   inj_p     = qs > 0;
   qs(inj_p) = qs(inj_p) .* comp(inj_p,1);
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_src(src)
   % Explicit sources defined by (e.g.) 'addSource'.

   qi = src.cell;
   qs = src.rate;

   % Injection sources have positive rate into reservoir.
   %
   in = find(src.rate > 0);
   if ~isempty(in),
      qs(in) = qs(in) .* src.sat(in,1);
   end
end

%--------------------------------------------------------------------------

function [qi, qs] = contrib_bc(G, resSol, bc, is_int)
   % Contributions from boundary conditions as defined by 'addBC'.

   qs = zeros([G.faces.num, 1]);
   dF = false([G.faces.num, 1]);

   isDir = strcmp('pressure', bc.type);
   isNeu = strcmp('flux',     bc.type);

   dF(bc.face(isDir))      = true;
   cfIx                    = dF(G.cellFaces(:,1));

   qs(G.cellFaces(cfIx,1)) = -resSol.cellFlux(cfIx);
   qs(bc.face(isNeu))      = bc.value(isNeu);

   % Injection BC's have positive rate (flux) into reservoir.
   %
   is_inj = qs > 0;
   if any(is_inj),
      qs(is_inj) = qs(is_inj) .* bc.sat(is_inj(bc.face), 1);
   end

   is_outer = ~is_int;

   qi = sum(G.faces.neighbors(is_outer,:), 2);
   qs = qs(is_outer);
end
