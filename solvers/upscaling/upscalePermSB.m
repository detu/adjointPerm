function perm = upscalePermSB(g, cg, rock, varargin)
% Compute upscaled Stokes-Brinkman permeabilites using flow based upscaling.
%
% SYNOPSIS:
%   perm_coarse = upscalePermSB(G, CG, rock, fluid)
%   perm_coarse = upscalePermSB(G, CG, rock, fluid, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid
%              model (CG) defined by function 'generateCoarseGrid'.
%
%   rock     - Rock data structure with valid field 'perm'.  The
%              permeability is assumed to be in measured in units of
%              metres squared (m^2).
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters. 
%
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
%   upscalePerm, extractSubgrid, computeMimeticIP.

%{
#COPYRIGHT#
%}

% $Id: upscalePermSB.m  $

   assert ( isfield(g, {'cartDims'}))

   opt = struct('Verbose', false);   
   opt = merge_options(opt, varargin{:});
   
   if opt.Verbose,
      h = waitbar(0, 'Computing upscaled Stokes-Brinkman permeabilities...');
      fprintf('Computing upscaled Stokes-Brinkman permeabilities... '), tic
   end
       
   fluid    = initSingleFluid('mu', 1000, 'rho', 1); % mu = 1
   
   dims     = find(g.cartDims);
   nBlocks  = cg.cells.num; 
   perm     = zeros(cg.cells.num, numel(g.cartDims));
   
   ups = @(cells, K_cells) upscale_block(extractSubgrid(g, cells), ...
                                         K_cells, dims, fluid);
   map.sub_cells = @(b) find(any(cg.cells.subCells(:,b), 2));
  
   
   % Compute upscaled permeability for each coarse block
   for b = 1 : nBlocks,
            
      cells = map.sub_cells(b);  % Fine-scale cells present in 'blk'
      K_cells = rock.perm(cells, :);
          
      perm(b,:) = ups(cells, K_cells); 
            
      if opt.Verbose, waitbar(b / nBlocks, h), end
   end
   if opt.Verbose,
      toc, close(h)
     
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function k = upscale_block(subG, K_cells, dims, fluid)
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
fluid.mu_eff = fluid.mu;
Dofs         = findCartDofs(subG);
rock.perm    = K_cells;
k            = zeros(1, numel(dims));

%{'LEFT', 'RIGHT', 'FRONT', 'BACK', 'BOTTOM', 'TOP'};
tags = [ 1 2; 3 4; 5 6]; 

bnd_f = any(subG.faces.neighbors==0,2);
ind   = bnd_f(subG.cellFaces(:,1));
faceAndTag = subG.cellFaces(ind, :);

for i = dims

   faces1 = faceAndTag(faceAndTag(:,2) == tags(i,1));
   faces2 = faceAndTag(faceAndTag(:,2) == tags(i,2));
   ind        = bnd_f;
   ind(faces1) = false;
   ind(faces2) = false;
   noflow_faces = find(ind);
   
   bc = addBCSB([], faces1, 'pressure', repmat(1*barsa, numel(faces1),  1), ...
                subG, Dofs);
   bc = addBCSB(bc, faces2, 'pressure', repmat(0, numel(faces2),  1), ...
                subG, Dofs);
   bc = addBCSB(bc, noflow_faces, 'velocity_n',...
                repmat(0, numel(noflow_faces),1), subG, Dofs);
   
   Ssb = makeSystemSB(subG, Dofs, rock, fluid, 'bc', bc);
   
   [Ssb, rSol] = solveSystemSB(Ssb, subG, Dofs, 'bc', bc);
   rSol        = nodeToCellData(rSol, subG, Dofs);
      
   area = sum(subG.faces.areas(faces2,:));
   L    = abs(subG.faces.centroids(faces1(1),i)...
           -subG.faces.centroids(faces2(1),i));

   q    = abs(sum(rSol.faceFlux(faces2)));  
   k(i) = q*L/(1*barsa*area);   
end
end

%--------------------------------------------------------------------------
