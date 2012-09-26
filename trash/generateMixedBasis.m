function CS = generateMixedBasis(G, S, CG, CS, weight, varargin) %#ok
% generateMixedBasis -- Construct coarse system interior basis functions or
%                       update basis functions for a list of faces. 
%
% SYNOPSIS:
%   basis = generateMixedBasis(G, S, CG, CS, weight)
%   basis = generateMixedBasis(G, S, CG, CS, weight, 'pn', pv, ...)
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
%   weight  - vector of size G.cells.num-by-1 giving weighting function
%             for generating the velocity basis functions.
%
%   'pn'/pv - List of 'key'/value pairs defining optional properties.  The
%             supported options are:
%               - Verbose -- Whether or not to emit progress reports while
%                            computing the flux basis functions.  Moreover,
%                            provide timing of the computations.
%                            Logical.  Default value = FALSE.
%
%               - Overlap -- Number of fine-grid cells in each physical
%                            direction with which to extend the supporting
%                            domain of any given basis function.
%                            Non-negative integer.  Default value = 0.
%
%               - updateFaces  -- list of faces to be updated. 
%                                 Default value = []. NB: must be inner
%                                 faces. 
%
%
% RETURNS:
%   basis   - Matrix of size nff-by-ncf
%               (nff = number of fine-scale faces,
%                ncf = number of coarse-scale faces)
%             of multiscale velocity basis functions.
%
%             To generate the coarse grid mass matrix, compute the product
%
%                    CS.B = basis.' * Do' * S.B * Do * basis
%
%             with S.B denoting the (up-to-date) mass matrix on the fine
%             grid, and Do is the matrix mapping face-fluxes to
%             half-face-fluxes.
%
% NOTE:
%   Using overlapping domains enables capturing more complex flow patterns,
%   particularly for very coarse grids, at the expense of increased
%   coupling in the resulting systems of linear equations.
%
%   Function generateMixedBasis is ONLY supported if the system structure
%   'S' was created by passing option pair ('Type','mixed') to function
%   assembleMimeticSystem·
%
% SEE ALSO:
%   generateBasis, mixedSymm, assembleMimeticSystem.

% $Id: generateMixedBasis.m 1734 2009-03-16 18:10:56Z bska $

if ~isfield(S, 'B'),
   error('generateMixedBasis:SystemType:NotSupported', ...
         ['Fine-grid inner product must be represented in ''mixed'' form.\n', ...
          'See option ''Type'' to function ''assembleMimeticSystem'''])
end

opt = struct('Verbose', false, 'Overlap', 0, 'updateFaces', []);
opt = merge_options(opt, varargin{:});

verbose = opt.Verbose;
overlap = opt.Overlap;
subCellsOverlap = CG.cells.subCells;
if overlap > 0,
   if verbose,
      fprintf('Generating overlap regions ... ')
      tic
   end
   % Generate neighbour matrix
   intF = prod(double(G.faces.neighbors), 2) > 0;
   neighborMatrix = sparse([(1 : G.cells.num) .';
                            double(G.faces.neighbors(intF, 1));
                            double(G.faces.neighbors(intF, 2))], ...
                           [(1 : G.cells.num) .';
                            double(G.faces.neighbors(intF, 2));
                            double(G.faces.neighbors(intF, 1))], ...
                           1, G.cells.num, G.cells.num);
   for j = 1 : overlap,
      subCellsOverlap = logical(neighborMatrix * subCellsOverlap);
   end
   tocif(verbose)
end

if verbose,
   fprintf('Computing flux basis functions ... ');
   tic
end

if ~isfield(CS, 'basis') 
   CS.basis = sparse(S.sizeD(2), CG.faces.num);
end

if isempty(opt.updateFaces)
   activeFaces    = find(prod(double(CG.faces.neighbors), 2)); 
else
   activeFaces = opt.updateFaces;   
end
   
numActiveFaces = numel(activeFaces);
numCF  = size(G.cellFaces,1);

cellNo = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';

% normal sign of each cell face.
sgn    = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo) - 1;

Do     = spdiags(sgn, 0, numCF, numCF)*S.D;

if verbose, h = waitbar(0, 'Computing flux basis functions ...'); end

for k = 1 : numActiveFaces,
   face   = activeFaces(k);
   blocks = CG.faces.neighbors(face,:);

   iG1 = CG.cells.subCells(:, blocks(1)); 
   iG2 = CG.cells.subCells(:, blocks(2));
   
   iG  = any(subCellsOverlap(:, blocks), 2);
   iF  = logical(S.C   * iG);
   iH  = logical(S.D.' * iF);
   
   % Compute synthetic source terms and define hybrid system.
   g0 = zeros([G.cells.num, 1]);
   sG1  = G.cells.volumes(iG1) .* weight(iG1);
   sG2  = G.cells.volumes(iG2) .* weight(iG2);

   g0(iG1) =   sG1 ./ sum(sG1);         % Source in block 1
   g0(iG2) = - sG2 ./ sum(sG2);         % Sink in block 2
   subG    = g0(iG);
   subF = zeros([nnz(iF), 1]);
   subH = zeros([nnz(iH), 1]);

   subB = S.B(iF,iF);
   subC = S.C(iF,iG);
   subD = S.D(iF,iH);
   
   
   neumannFaces = (sum(subD) == 1) .';
  
   
   subDo = Do(iF, iH);

% Solve (symmetric) hybrid system, no-flow BC's.
% remove all non-neumann faces from D and H
   fv = mixedSymm(subB, subC, subD(:,neumannFaces), subF, subG, ...
                  subH(neumannFaces), subDo, 'Regularize', true);
   
   CS.basis(:, face) = sparse(find(iH), 1, fv, S.sizeD(2), 1);
   if verbose, waitbar(k / numActiveFaces, h); end
end

if verbose, toc, close(h); end
