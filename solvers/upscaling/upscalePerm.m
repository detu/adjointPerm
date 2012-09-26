function perm = upscalePerm(g, cg, rock, varargin)
% Compute upscaled permeabilites using flow based upscaling.
%
% SYNOPSIS:
%   perm_coarse = upscalePerm(G, CG, rock, fluid)
%   perm_coarse = upscalePerm(G, CG, rock, fluid, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid
%              model (CG) defined by function 'generateCoarseGrid'.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%              permeability is assumed to be in measured in units of
%              metres squared (m^2).
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters. 
%
%                - S        -- Mimetic linear system data structures as 
%                              defined by function 'computeMimeticIP'.
%
%                - LinSolve --
%                         Handle to linear system solver software to which
%                         the fully assembled system of linear equations
%                         will be passed.  Assumed to support the syntax
%
%                             x = LinSolve(A, b)
%
%                         in order to solve a system Ax=b of linear eq.
%                         Default value: LinSolve = @mldivide (backslash).
%                - Verbose -- 
%                         Whether or not to emit progress reports during
%                         the assembly process.
%                         Logical.  Default value = FALSE.
%
% RETURNS:
%   perm_coarse - Upscaled permeability field
%
% REMARK:
%   The flow based upscaling applied here is based on the assumption that
%   the grid is Cartesian. 
%
% SEE ALSO:
%   extractSubgrid, computeMimeticIP.

%{
#COPYRIGHT#
%}

% $Id: upscalePerm.m  $

   assert ( isfield(g, {'cartDims'}))

   opt = struct('Verbose', false, 'LinSolve',   @mldivide, 'S', []);   
   opt = merge_options(opt, varargin{:});
   
   if opt.Verbose,
      h = waitbar(0, 'Computing upscaled permeabilities...');
      fprintf('Computing upscaled permeabilities... '), tic
   end
   
   % make mimetic ip if not supplied
   if isempty(opt.S)   
      s = computeMimeticIP(g, rock);
   else
      s = opt.S;
   end
      
   fluid    = initSingleFluid('mu', 1000, 'rho', 1); % mu = 1
   
   dims     = find(g.cartDims);
   nBlocks  = cg.cells.num; 
   perm     = zeros(cg.cells.num, numel(g.cartDims));
   
   ups = @(cells, subS) upscale_block(extractSubgrid(g, cells), subS, ...
                                      dims,  fluid, opt.LinSolve);
   
   map.sub_cells = @(b) find(any(cg.cells.subCells(:,b), 2));
   hfix          = cumsum([0; double(g.cells.numFaces)]);
   map.hf        = @(c) mcolon(hfix(c) + 1, hfix(c + 1)).';

   subS.ip   = s.ip;
   subS.type = s.type;
   
   % Compute upscaled permeability for each coarse block
   for b = 1 : nBlocks,
            
      cells = map.sub_cells(b);  % Fine-scale cells present in 'blk'
      hf    = map.hf(cells);     % Fine-scale half-faces in 'blk'
      
      subS.BI = s.BI(hf,hf);     % Extract subsystem 
        
      perm(b,:) = ups(cells, subS); 
            
      if opt.Verbose, waitbar(b / nBlocks, h), end
   end
   if opt.Verbose,
      toc, close(h)
     
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function k = upscale_block(subG, subS, dims, fluid, linsolve)
%  
%  Compute perm in x-direction:
%
%                  v*n = 0    
%              ______________  
% l_boundary  |              | r_boundary          
%             |              |
%     p1 = 1  |              | p2 = 0
%             |______________|
%                   v*n = 0
% 
%             |--------------|    
%                    L  
%                    
%           grad p = (p2 - p1)/L 
%           q1: flux over right boundary
%           kx* = (q1*L)/(area(r_boundary)*(p2 - p1))
%     

k = zeros(1, numel(dims));
rSol = initResSol(subG, 0.0);

%{'LEFT', 'RIGHT', 'FRONT', 'BACK', 'BOTTOM', 'TOP'};
tags = [ 1 2; 3 4; 5 6]; 

bnd_f = any(subG.faces.neighbors==0,2);
ind = bnd_f(subG.cellFaces(:,1));
faceAndTag = subG.cellFaces(ind, :);

for i = dims

   faces1 = faceAndTag(faceAndTag(:,2) == tags(i,1));
   faces2 = faceAndTag(faceAndTag(:,2) == tags(i,2));

   bc = addBC([], faces1, 'pressure', 1*barsa);
   bc = addBC(bc, faces2, 'pressure', 0);

   rSol = solveIncompFlow(rSol, [], subG, subS, fluid, 'bc', bc, ...
                          'LinSolve',linsolve);
   area = sum(subG.faces.areas(faces2,:));
   L    = abs(subG.faces.centroids(faces1(1),i)...
           -subG.faces.centroids(faces2(1),i));

   q    = abs(sum(rSol.faceFlux(faces2)));  
   k(i) = q*L/(1*barsa*area);  
end
end

%--------------------------------------------------------------------------
