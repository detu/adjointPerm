function [V, P] = evalBasisFuncMixed(faces, g, cg, B, C, D, w, mob, varargin)
% evalBasisFuncMixed -- Compute multiscale mixed basis functions for
%                       selected faces.
%
% SYNOPSIS:
%   [V, P] = evalBasisFuncMixed(faces, G, CG, B, C, D, weight, Lt)
%   [V, P] = evalBasisFuncMixed(faces, G, CG, B, C, D, weight, Lt, ...
%                               'pn1', pv1, ...)
%
% PARAMETERS:
%   faces    - List (array) of coarse faces for which (new) basis function
%              values are requested.  Must be a list of explicit indices.
%              Function 'evalBasisFuncMixed' does not support a LOGICAL
%              array.
%
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid model
%              (CG) defined by function 'generateCoarseGrid'.
%
%   B, C, D  - Component system matrices from underlying fine-grid
%              (mimetic) mixed hybrid discretisation of incompressible
%              pressure equation. 
%
%   w        - Array of numerically evaluated synthetic volume source
%              weighting term.  One scalar value for each cell in the
%              underlying fine grid g.  All values must be supplied, even
%              if some of the cells do not participate in the any of the
%              basis functions for 'faces'.
%
%   mob      - Total mobility.  One scalar value for each cell in the
%              underlying fine model.
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
%                - bc  -- External boundary conditions.
%                         Must be a boundary condition data structure as
%                         defined by function 'addBC'.  Default value is []
%                         (an empty array), meaning that no boundary
%                         conditions are present in the model.
%
% RETURNS:
%   V - Matrix of size nfhf-by-nchf of multiscale flux basis functions.
%       Here 'nfhf' is the number of fine-scale half-faces (half-contacts)
%       and 'nchf' is the number of coarse-scale half-faces in the set
%       represented by 'faces'.  Note that the basis function values are
%       really 'B' times the actual flux values with 'B' denoting the
%       fine-scale mass matrix.
%
%   P - Matrix of size nfc-by-nchf of pressure basis functions (computed
%       using a fluid with total mobility equal to one).  Here, 'nfc' is
%       the number of cells in the fine-scale model.
%
% EXAMPLE:
%   % Define simple geometry and rock properties.
%   G = computeGeometry(cartGrid([50, 50, 4]));
%   rock = struct('perm' , repmat(100, [g.cells.num, 1]), ...
%                 'poros', repmat(0.3, [g.cells.num, 1]));
%
%   % Define coarse grid.
%   CG = generateCoarseGrid(G, partitionUI(G, [5, 5, 1]));
%
%   % Assemble hybrid pressure system matrices.
%   S = assembleMimeticSystem(G, rock, 'type', 'mixed');
%
%   % Compute basis functions for all inner coarse faces using cell
%   % porosity as a weighting function.
%   %
%   faces  = find(all(g.faces.neighbors > 0, 2));
%   mob    = ones([G.cells.num, 1])
%   [V, P] = evalBasisFuncMixed(faces, G, CG, S.B, S.C, S.D, ...
%                               rock.poro, mob );
%
% SEE ALSO:
%   generateCoarseGrid, generateCoarseSystem, assembleMimeticSystem.
%
% NOTE:
%   Does currently not support overlap [TODO].

% $Id: evalBasisFuncMixed.m 1829 2009-03-23 14:01:07Z bska $

   if ~isnumeric(faces),
      error(id('Faces:NonNumeric'), ...
            'Argument ''faces'' must be a numeric array.');
   end

   opt = struct('Verbose', false, 'src', [], 'bc', []);
   opt = merge_options(opt, varargin{:});

   % Build various index mappings.
   if ~isempty(opt.bc),
      [nsub, sub] = subFaces(g, cg);
      sub_ix      = cumsum([0; nsub]);
   end

   cellno  = rldecode((1 : g.cells.num).', double(g.cells.numFaces));
   % normal sign of each cell face.
   sgn   = 2*double(g.faces.neighbors(g.cellFaces(:,1), 1) == cellno) - 1;

   part    = get_partition(cg);
   ncells  = accumarray(part, 1);
   % number of half faces in a block 
   nhf_blk = accumarray(part, double(g.cells.numFaces));
   % ~= nsub!
   hfix    = cumsum([0; double(g.cells.numFaces)]);
   ffix    = sparse(double(g.cellFaces(:,1)), part(cellno), 1, ...
                    g.faces.num, cg.cells.num);

   % Driving forces (sources and BCs).
   [f_bc, h_bc, is_dirichlet] = expand_bc(g, opt.bc);
   theta = get_coarse_weighting(g, part, w, opt.src);

   % Initialize matrices used to construct basis (i-index, j-index, value)
   iv = []; jv = []; sv = [];
   ip = []; jp = []; sp = [];
   

   numCF = size(g.cellFaces,1);
   Do    = spdiags(sgn, 0, numCF, numCF)*D;

   if opt.Verbose,
      k = 0;
      h = waitbar(k, 'Computing flux and pressure basis functions...');
      nFaces = numel(faces);
      fprintf('Computing flux and pressure basis functions... '), tic
   end
   
   % Make into row vector for benefit of FOR loop.
   faces = reshape(faces, 1, []);
   for face = faces,
      blk  = cg.faces.neighbors(face,:); % Blocks connected to 'face'.
      blk  = blk(blk > 0);
      nblk = numel(blk);                 % Number of blocks sharing 'face'.
      nc   = ncells(blk);                 
      % if blk = [blk1, blk2], mod makes sure that the indices are correct
      % for blk2 also
      iG = mod(find(cg.cells.subCells(:, blk)) - 1, g.cells.num) + 1;           
      iF = mcolon(hfix(iG) + 1, hfix(iG + 1)) .'; % mcolon yields row vector.
      iH = find(any(ffix(:,blk), 2));
                      
      sB = B(iF,     iF );  sC = C(iF, iG); sD = D(iF,      iH );
      sF = zeros(size(iF)); sG = theta(iG); sH = zeros(size(iH));

      % Update for effects of several phases...
      
      dmob = spdiags(1./mob(cellno(iF)), 0, numel(iF), numel(iF));
      sB  = dmob * sB;

      is_dir = false(size(iH));
      lam    = zeros(size(iH));

      if nblk > 1,
         % Put sink in block 2.
         %
         sG(nc(1) + 1 : end) = -sG(nc(1) + 1 : end);       
      else
         % The coarse face 'face' is on the boundary of the domain.
         % Consequently, we need to include boundary conditions for 'face'
         % but not for any other external coarse faces which might happen
         % to connect to 'blk'.  Such other faces always get no-flow
         % conditions (homogeneous Neumann).
         %
         ih      = sub(mcolon(sub_ix(face) + 1, sub_ix(face + 1)));
         is_bdry = false(size(iH));
         is_bdry(ismember(iH, ih)) = true;

         is_dir(is_bdry) = is_dirichlet(iH(is_bdry));
         is_neu          = is_bdry & ~is_dir;

         % This code assumes that a coarse face is either Dirichlet or
         % Neumann (not both).  That assumption must be revisited if the
         % following assertion fails.  We also fail if the face is neither
         % Dirichlet nor Neumann.
         %
         assert (xor(sum(is_dir) > 0, sum(is_neu) > 0));

         %lam(is_dir) = f_bc(iH(is_dir));
         sF          = sF - sD(:,is_dir)*lam(is_dir);
         sH(is_neu)  = h_bc(iH(is_neu));
         denom       = sum(sH(is_neu));
         if abs(denom) > sqrt(eps(denom)),
            sH(is_neu) = sH(is_neu) / denom;
         end
      end

      do_reg = sum(is_dir) == 0;  % Need to set pressure zero level?

      sDo = Do(iF, iH);
      
      % identify neumann faces
      neuFaces = (sum(sD) == 1) .';

      [fv, p] = mixedSymm(sB, sC, sD(:,neuFaces), ...
                          sF, sG, sH(  neuFaces), ...
                          sDo, 'Regularize', do_reg);
      
      % convert from faceFlux to cellFlux and multiply by dmob            
      v  = sD * fv .* sgn(iF);        
      v  = dmob * v; 
      j1 = zeros(size(v)); j1(1) = 1;
      j2 = zeros(size(p)); j2(1) = 1;
      if nblk > 1,
         v(nhf_blk(blk(1)) + 1 : end) = -v(nhf_blk(blk(1)) + 1 : end);
         p(         nc(1)  + 1 : end) = -p(         nc(1)  + 1 : end);

         j1(nhf_blk(blk(1)) + 1)      = 1;
          
         j2(nc(1)           + 1)      = 1;

         a = [p(        1 : nc(1)).' * sG(        1 : nc(1)); ...
              p(nc(1) + 1 :  end ).' * sG(nc(1) + 1 :  end )];
         a = rldecode([a(1); -a(2)], nc(:));
      else
         a = repmat(p.' * sG, [nc, 1]);
      end

      p = p - a;  % Orthogonalize against source term 'sG' (per coarse blk)

      

      if abs(p' * sG) > 2 * numel(p) * eps(norm(p,inf)),
         warning(msgid('Orthogonality:Questionable'),       ...
                 ['Questionable orthogonality of p-basis ', ...
                  'for face ''%d'''], face);
      end
      iv = [iv; iF]; jv = [jv; j1]; sv = [sv; v]; %#ok
      ip = [ip; iG]; jp = [jp; j2]; sp = [sp; p]; %#ok
      

      if opt.Verbose, k = k + 1; waitbar(k / nFaces, h), end
   end

   if opt.Verbose, toc, fprintf('Assembling output matrices ... '), tic, end

   % Assemble final output.
   jv = cumsum(jv);            % => max(jv) == jv(end)
   jp = cumsum(jp);            % => max(jp) == jp(end)  
   V = sparse(iv, jv, sv, size(g.cellFaces,1), jv(end));
   P = sparse(ip, jp, sp, g.cells.num        , jp(end));
   if opt.Verbose, toc, close(h), end
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

%--------------------------------------------------------------------------

function [f, h, d] = expand_bc(g, bc)
   f = zeros([g.faces.num, 1]);
   h = zeros([g.faces.num, 1]);
   d = false([g.faces.num, 1]);

   if ~isempty(bc),
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      is_dir = strcmp('pressure', bc.type);
      f(bc.face(is_dir)) = bc.value(is_dir);
      d(bc.face(is_dir)) = true;

      is_neu = strcmp('flux', bc.type);
      h(bc.face(is_neu)) = - bc.value(is_neu);
   end
end

%--------------------------------------------------------------------------

function s = id(s)
   s = ['evalBasisFuncMixed:', s];
end
