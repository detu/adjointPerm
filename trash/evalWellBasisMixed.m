function [V, P, Q] = evalWellBasisMixed(well, g, cg, B, C, D, w, varargin)
% evalWellBasis -- Compute multiscale basis functions for single well.
%
% SYNOPSIS:
%   [V, P, Q] = evalWellBasisMixed(w, G, CG, B, C, D, weight)
%   [V, P, Q] = evalWellBasisMixed(w, G, CG, B, C, D, weight, 'pn1', pv1, ...)
%
% PARAMETERS:
%   w        - Single well structure as defined by functions 'addWell' and
%              'assembleWellSystem'.
%
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid model
%              (CG) defined by function 'generateCoarseGrid'.
%
%   B, C, D -  Component system matrices from underlying fine-grid
%              (mimetic) mixed hybrid discretisation of incompressible
%              pressure equation.  
%
%   weight   - Array of numerically evaluated synthetic volume source
%              weighting term.  One scalar value for each cell in the
%              underlying fine grid G.  All values must be supplied, even
%              if some of the cells do not participate in the any of the
%              basis functions for 'faces'.
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.  The
%              supported options are:
%                - src -- Explicit source terms in the underlying fine grid
%                         model which must be taken into account when
%                         generating the basis functions.  Note that basis
%                         functions in blocks containing an explicit source
%                         term will be generated based solely on the
%                         explicit source.  The values in 'weight'
%                         pertaining to such blocks will be ignored.
%
%                         Must be a source data structure as defined by
%                         function 'addSource'.  Default value is [] (an
%                         empty array), meaning that no explicit sources
%                         are present in the model.
%
% RETURNS:
%   V - Matrix of size nfhf-by-ncb of multiscale well-to-block flux basis
%       functions.  Here 'nfhf' is the number of fine-scale half-faces
%       (half-contacts) and 'ncb' is the number of coarse blocks
%       intersected by the well 'w'.  Note that the basis function values
%       are the actual flux values (as oposed to in the hybrid case).
%
%   P - Matrix of size nfc-by-ncb of (well-to-block) pressure basis
%       functions (computed using a fluid with total mobility equal to
%       one).  Here, 'nfc' is the number of cells in the fine-scale model.
%
%   Q - Matrix of size nwc-by-ncb of (well-to-block) well rate basis
%       functions.  Here, 'nwc' is the number of perforations (i.e.,
%       cells intersected by) the well 'w'.
%
% SEE ALSO:
%   evalBasisFunc, generateCoarseWellSystem.

% $Id: evalWellBasisMixed.m 1820 2009-03-23 12:41:49Z bska $

   opt = struct('Verbose', false, 'src', []);
   opt = merge_options(opt, varargin{:});

   cellno = rldecode((1 : g.cells.num).', double(g.cells.numFaces));
   % normal sign of each cell face.
   sgn   = 2*double(g.faces.neighbors(g.cellFaces(:,1), 1) == cellno) - 1;

   part   = get_partition(cg);
   hfix   = cumsum([0; double(g.cells.numFaces)]);
   ffix   = sparse(double(g.cellFaces(:,1)), part(cellno), 1, ...
                   g.faces.num, cg.cells.num);

   theta = get_coarse_weighting(g, part, w, opt.src);

   iv = []; jv = []; sv = [];
   ip = []; jp = []; sp = [];
   iq = []; jq = []; sq = [];

   pwc = part(well.cells);
   blk = reshape(unique(pwc), 1, []);  % Row vector needed in FOR loop.
   for b = blk,
      iG = mod(find(cg.cells.subCells(:, b)) - 1, g.cells.num) + 1;
      iF = mcolon(hfix(iG) + 1, hfix(iG + 1)) .'; % mcolon -> row vector.
      iH = find(any(ffix(:,b), 2));

      % Reservoir discrete system from block 'b'.
      sB = B(iF,     iF );  sC =  C(iF, iG); sD = D(iF,      iH );
      sF = zeros(size(iF)); sG = -theta(iG); sH = zeros(size(iH));

      % Append well block 'b' well system to reservoir discrete system.
      iFw = find(pwc == b);

      sB  = blkdiag(sB,  well.S.B(iFw,iFw));
      sC  = vertcat(sC , well.S.C (iFw,iG ));
      % force p_wf = 0
      DHat = blkdiag(sD, diag(ones( nnz(iFw), 1)));
      sD  = blkdiag(sD , sparse(numel(iFw), 0));
      sF  = vertcat(sF , zeros([numel(iFw), 1]));

      orientation = [sgn(iF); ones(numel(iFw), 1)];
      n  = length(orientation);   
      % blir denne riktig størrelse?
      Do =  spdiags(orientation, 0, n, n)*DHat;
      
      % identify neumann faces
      neuFaces = (sum(sD) == 1) .';
      [v, p] = mixedSymm(sB, sC, sD(:, neuFaces), sF, sG, sH(neuFaces), Do);
                 
      ixr = 1 : numel(iF);
      ixw = numel(iH) + 1 : numel(v);
      
      %convert reservoir fluxes from faceFlux to cellFlux
      f = sD(ixr, :)*v(1:numel(iH)).*sgn(iF);  
     
      p = p + (p' * sG);         % Orthogonalize aginst src.      
     
      if abs(p' * sG) > 2 * numel(p) * eps(norm(p,inf)),
         warning(msgid('Orthogonality:Questionable'),       ...
                 ['Questionable orthogonality of p-basis ', ...
                  'for face ''%d'''], face);
      end
      % Package solutions (i.e. the well basis functions) into sparse
      % matrix format.
      j1 = zeros(size(iF ));  j1(1) = 1;
      j2 = zeros(size(iG ));  j2(1) = 1;
      j3 = zeros(size(iFw));  j3(1) = 1;

      iv = [iv; iF ]; jv = [jv; j1]; sv = [sv; f]; %#ok
      ip = [ip; iG ]; jp = [jp; j2]; sp = [sp; p];       %#ok
      iq = [iq; iFw]; jq = [jq; j3]; sq = [sq; -v(ixw)]; %#ok
   end

   % Assemble final solution output.  Note that, intially, ALL(j? >= 0) so
   % following the call to CUMSUM, MAX(j?) == j?(end).
   %
   jv = cumsum(jv);  V = sparse(iv, jv, sv, size(B,1), jv(end));
   jp = cumsum(jp);  P = sparse(ip, jp, sp, size(C ,2), jp(end));
   jq = cumsum(jq);  Q = sparse(iq, jq, sq, numel(pwc), jq(end));
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function p = get_partition(cg)
   p     = zeros([size(cg.cells.subCells,1), 1]);
   [i,j] = find(cg.cells.subCells);
   p(i)  = j;
end

%--------------------------------------------------------------------------

function theta = get_coarse_weighting(g, p, w, src)
   % Initially, assume there are no explicit/external sources.
   %
   theta = w .* g.cells.volumes;

   if ~isempty(src),
      % Update for explicit sources if there nevertheless are some...

      % Determine coarse blocks already containing external sources.
      %
      has_src = unique(p(src.cell));

      % Eliminate previous (synthetic) weighting and apply correct source.
      %
      theta(ismember(p, has_src)) = 0;
      theta(src.cell)             = src.rate;
   end

   % Note:
   %   We need to normalize the (synthetic or explicit) source term 'theta'
   %   such that \int_{B_i} theta d\Omega == 1 lest the basis functions be
   %   inconsistent.
   %
   denom = accumarray(p, theta);
   theta = theta ./ denom(p);
end
