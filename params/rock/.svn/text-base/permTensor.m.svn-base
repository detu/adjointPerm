function [K, r, c] = permTensor(rock, dim)
%Expand permeability tensor to full format.
%
% SYNOPSIS:
%    K        = permTensor(rock, dim)
%   [K, r, c] = permTensor(rock, dim)
%
% PARAMETERS:
%   rock - Rock data structure with valid field 'perm'--one value of the
%          permeability tensor for each cell (nc) in the discretised model.
%
%          The field 'rock.perm' may have ONE column for a scalar
%          (isotropic) permeability in each cell, TWO or THREE columns for
%          a diagonal permeability in each cell (in two or three space
%          dimensions, respectively) and THREE or SIX columns for a
%          symmetric, full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%   dim  - Number of space dimensions in the discretised model.  Must be
%          two or three.
%
% RETURNS:
%   K    - Expanded permeability tensor.  An nc-by-(dim^2) array of
%          permeability values as described above.
%
%   r, c - Row- and column indices from which the two-form
%
%                (x, Ky) = \sum_r \sum_c x(r) * K(r,c) * y(c)
%
%          may be easily evaluated by a single call to SUM (see EXAMPLE).
%          OPTIONAL.  Only returned if specifically requested.
%
% EXAMPLE:
%   % Compute n'Kg gravity term on each cell face (half contact).
%   %
%   cellNo = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2).';
%   sgn    = 2*double(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo) - 1;
%   % cfn  = cell-face normal = *outward* face normals on each cell.
%   cfn    = bsxfun(@times, G.faces.normals(G.cellFaces(:,1), :), sgn);
%
%   dim    = size(G.nodes.coords, 2);
%   [K, r, c] = permTensor(rock, dim);
%
%   g   = gravity();   g = reshape(g(1 : dim), 1, []);
%   nKg = sum(cfn(:,r) .* bsxfun(@times, K(cellNo,:), g(c)), 2);
%
% SEE ALSO:
%   computeMimeticIP, computeTrans, solveBlackOilWellSystem.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-08 13:06:46 +0200 (ma, 08 jun 2009) $
% $Revision: 2348 $

   if isempty(rock) || ~isfield(rock, 'perm'),
      error(msgid('Rock:Empty'), ...
            '''rock'' structure must contain valid field ''perm''.')
   end

   [nc, nk] = size(rock.perm);

   if dim == 2,
      switch nk,
         case 1,
            % Isotropic.

            K = rock.perm * [1, 0, 1];

         case 2,
            % Diagonal.

            K = [rock.perm(:,1), zeros([nc,1]), rock.perm(:,2)];

         case 3,
            % Full, symmetric, 2-by-2.

            K = rock.perm;

         otherwise,
            error(msgid('PermDim:Unsupported'), ...
                 ['%d-component permeability value is not supported ', ...
                  'in two space dimensions.'], nk);
      end

      K = K(:, [1, 2, 2, 3]);
      r =      [1, 1, 2, 2] ;
      c =      [1, 2, 1, 2] ;

   elseif dim == 3,
      switch nk,
         case 1,
            % Isotropic.

            K = rock.perm * [1, 0, 0, 1, 0, 1];

         case 3,
            % Diagonal.

            K = [rock.perm(:,1), zeros([nc,2]), ...
                 rock.perm(:,2), zeros([nc,1]), ...
                 rock.perm(:,3)];

         case 6,
            % Full, symmetric, 3-by-3.

            K = rock.perm;

         otherwise,
            error(msgid('PermDim:Unsupported'), ...
                 ['%d-component permeability value is not supported ', ...
                  'in three space dimensions.'], nk);
      end

      K = K(:, [1, 2, 3, 2, 4, 5, 3, 5, 6]);
      r =      [1, 1, 1, 2, 2, 2, 3, 3, 3] ;
      c =      [1, 2, 3, 1, 2, 3, 1, 2, 3] ;

   else
      error(msgid('PermDim:Unsupported'), ...
            '%d space dimensions are unsupported at this time.', dim);
   end
end
