function CS = generateBasis(G, S, CG, CS, weight, varargin)
% generateBasis -- Construct/update coarse system interior basis functions
%
% SYNOPSIS:
%   CS = generateBasis(G, S, CG, CS, blocks, weight)
%   CS = generateBasis(G, S, CG, CS, blocks, weight, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure.
%
%   S       - System struture.
%
%   CG      - Coarse grid structure.
%
%   CS      - Coarse system structure.
%
%   weight  - Vector giving weighting function for generating 
%             the velocity basis functions.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Verbose -- Whether or not to emit progress reports while
%                            computing basis functions as well as timing
%                            the computations.
%                            Logical.  Default value = FALSE.
%
%               - ActiveBndFaces -- Vector of boundary faces that are
%                                   active (not no-flow) and are to be
%                                   assigned a basis function.
%                                   Default value = [].
%
%               - Blocks -- Coarse blocks for which basis functions will be
%                           generated.  Can be given as [] if only basis
%                           for ActiveBndFaces are to be generated.
%                           Default value = 1:CG.cells.num (generate basis
%                           for all blocks).
%
% RETURNS:
%   CS - System structure with added/updated field:
%       -basis    - Matrix of size nfhf-by-ncf
%                   (nfhf = number of fine-scale half-faces,
%                    ncf  = number of coarse-scale faces)
%                   of multiscale velocity basis functions.  Note that the
%                   basis function values are really `B' times the actual
%                   velocity values with `B' denoting the fine-scale mass
%                   matrix.
%
%             To generate the coarse grid mass matrix, compute the
%             product
%
%                    CS.B = CS.basis.' * S.BI * CS.basis
%
%             with S.BI denoting the (up-to-date) inverse mass matrix on
%             the fine grid (i.e. S.BI = INV(B)).
%
% NOTE:
%   Function generateBasis is NOT supported if the system structure 'S'
%   was created by passing option pair ('Type','mixed') to function
%   assembleMimeticSystem.
%
% SEE ALSO:
%   generateMixedBasis, schurComplementSymm, assembleMimeticSystem,
%   putBCMS.

% $Id: generateBasis.m 1705 2009-03-12 16:23:43Z bska $

if ~isfield(S, 'BI'),
   error(id('SystemType:NotSupported'), ...
         ['Fine-grid inner product must be represented in ''hybrid'' form.\n', ...
          'See option ''Type'' to function ''assembleMimeticSystem'''])
end
opt = struct('Verbose', false,         ...
             'Blocks', 1:CG.cells.num, ...
             'ActiveBndFaces', []);
opt = merge_options(opt, varargin{:});

verbose   = opt.Verbose;
activeBnd = opt.ActiveBndFaces;
blocks    = opt.Blocks;

% Check activeBnd is on boundary (needed if function is used to update
% basis)
if sum(prod(double(CG.faces.neighbors(activeBnd)),2)) ~= 0,
    error(id('ActiveBndFaces:NotBoundary'), ...
         'Given ''activeBnd'' are not on boundary.');
end

if verbose,
   fprintf('Computing flux basis functions... ');
   tic
end

% Find coarse faces of the blocks given in the input parameter "blocks"
blockFaces = unique(CG.cellFaces(ismember(CG.cellFaces(:,1),blocks),2));

% Check if CS.basis is initialized and do so if not
if ~(isfield(CS,{'basis'}) && issparse(CS.basis) && ...
   all(size(CS.basis)==[S.sizeB(1) length(CG.cellFaces)]))
   CS.basis = sparse(S.sizeB(1), length(CG.cellFaces));
end

BIndex = S.C  * CG.cells.subCells;
DIndex = S.D' * BIndex;

% Extract faces to be updated
innerBlockFaces = blockFaces(prod(double(CG.faces.neighbors(blockFaces)),2)>0);
activeFaces = [innerBlockFaces; activeBnd];  
nFaces = length(activeFaces);

if verbose, h = waitbar(0, 'Computing flux basis functions...'); end

for k = 1 : nFaces,
   face  = activeFaces(k);
   cells = CG.faces.neighbors(face,:);
   boundary = prod(double(cells),2) == 0;

   if ~boundary,
      [iF1, iF2, iF] = getIx(BIndex(:, cells));
      [iG1, iG2, iG] = getIx(CG.cells.subCells(:, cells));
      [iH,  iH,  iH] = getIx(DIndex(:, cells));  % No typo!

      % Compute source terms and define hybrid system.
      % If a block does not contain a source, a syntetic source term is
      % used.
      g0   = zeros([G.cells.num, 1]);
      if nnz(S.RHS.g_src(iG1)) == 0,       % -> Synthetic source
         sG1 = weight(iG1);
      else
         sG1 = S.RHS.g_src(iG1);
      end
      if nnz(S.RHS.g_src(iG2)) == 0,       % -> Synthetic source
         sG2 = weight(iG2);
      else
         sG2 = S.RHS.g_src(iG2);
      end
      sG1 = sG1 .* G.cells.volumes(iG1);
      sG2 = sG2 .* G.cells.volumes(iG2);

      g0(iG1) =   sG1 ./ sum(sG1);         % Source in block 1
      g0(iG2) = - sG2 ./ sum(sG2);         % Sink in block 2

      subF = zeros(size(iF));
      subH = zeros(size(iH));
      subD = S.D  (iF,iH);

      dirF = false(size(subH));            % Empty if Neumann BC.
      dirCond = [];
   else %boundary face
      cell = nonzeros(double(cells));
      iF   = find(BIndex(:, cell));
      iG   = find(CG.cells.subCells(:, cell));
      iH   = find(DIndex(:, cell));

      % Compute source terms and define hybrid system.
      g0 = zeros([G.cells.num, 1]);

      if nnz(S.RHS.g_src(iG)) == 0,     % -> Synthetic source
         sG = weight(iG);
      else
         sG = S.RHS.g_src(iG);
      end

      sG     = sG .* G.cells.volumes(iG);
      g0(iG) = sG ./ sum(sG);           % Source in block

      subF = zeros(size(iF));
      dirF = logical(sparse(numel(iH), 1));
      subH = zeros(size(iH));
      subD = S.D(iF,iH);
      
      % Find faces on boundary face
      fineFInx  = iH(G.faces.tag(iH) == CG.faces.tag(face));
      localFInx = ismember(iH, fineFInx);

      % Extract info from fine system S
      if sum(S.RHS.h_bc(fineFInx))~=0             
         subH(localFInx)  =  S.RHS.h_bc(fineFInx)/sum(S.RHS.h_bc(fineFInx));             
      end
      dirF(localFInx)  =  S.RHS.dirichletFaces(fineFInx);
      dirCond          = -S.RHS.f_bc(iF(logical(subD * dirF)));

      localCFInx       = logical(sum(subD(:,localFInx),2));
      fineCFInx        = iF(localCFInx);
      subF(localCFInx) = S.RHS.f_bc(fineCFInx);
   end
   
   subG  = g0(iG);
   subBI = S.BI(iF,iF);
   subC  = S.C (iF,iG);

   % Solve (symmetric) hybrid system.
   % No-flow BC's for blocks away from boundary,
   % otherwise use fine system BC.
   lam       = zeros(size(iH));
   lam(dirF) = dirCond;
   do_reg    = sum(dirF) == 0;  % Need to set pressure zero level?
   [v, p, lam(~dirF)] = schurComplementSymm(subBI, subC, subD(:,~dirF), ...
                                            subF , subG, subH(  ~dirF), ...
                                            'regularize', do_reg);

   Bv = sparse(iF, 1, subC*p - subD*lam, S.sizeB(1), 1);

   cfIx = find(CS.D(:, face));
   if ~boundary, %inner face
      % 'left' and 'right' basis functions respectively
      %
      sgn = 2*(cells(1) == CG.faces.neighbors(face,1)) - 1;
      CS.basis(:, cfIx(1)) = sparse(iF1, 1,  sgn .* Bv(iF1), S.sizeB(1), 1);
      CS.basis(:, cfIx(2)) = sparse(iF2, 1, -sgn .* Bv(iF2), S.sizeB(1), 1);

   else %boundary
      CS.basis(:, cfIx(1)) = sparse(iF, 1,  Bv(iF), S.sizeB(1), 1);
   end
   if verbose, waitbar(k / nFaces, h); end
end

if verbose, toc, close(h); end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [ix1, ix2, ix] = getIx(v)
ix1 = find(v(:,1));
ix2 = find(v(:,2));
ix  = find(sum(v,2));
