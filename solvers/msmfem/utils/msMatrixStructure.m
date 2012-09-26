function A = msMatrixStructure(G, CG, varargin)
%Build synthetic matrix with same sparsity as hybrid coarse system matrix.
%
% SYNOPSIS:
%   A = msMatrixStructure(G, CG)
%   A = msMatrixStructure(G, CG, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure as described by 'grid_structure'.
%   CG      - Coarse grid data structure as defined by function
%             'generateCoarseGrid'.
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - bc -- External boundary conditions.
%                       Must be a boundary condition data structure as
%                       defined by function 'addBC'.  Default value is []
%                       (an empty array), meaning that no boundary
%                       conditions are present in the model.
%
% RETURNS:
%   A - A synthetic matrix with the same sparsity structure as the coarse
%       hybrid system matrix derived in solver function
%       'solveIncompFlowMS'.  This matrix is mostly usable for passing
%       directly to sparsity inspection function 'spy'.  In particular, the
%       matrix 'A' does not have any of the analytic properties (e.g.,
%       eigenvalues) of the coarse hybrid system matrix.
%
% LIMITATION:
%   This function does not at present support well connections.
%
% SEE ALSO:
%   generateCoarseSystem, solveIncompFlowMS, spy.

%{
#COPYRIGHT#
%}

% $Id: msMatrixStructure.m 1953 2009-03-31 10:54:12Z bska $

   opt = struct('bc', []);
   opt = merge_options(opt, varargin{:});

   % Determine, and enumerate, the active coarse (reservoir) faces.
   is_act = false([CG.faces.num, 1]);
   is_act(all(CG.faces.neighbors > 0, 2))    = true;
   is_act(bdry_supports_flow(G, CG, opt.bc)) = true;

   renum         = zeros([CG.faces.num, 1]);
   renum(is_act) = 1 : sum(double(is_act));

   % Extract active coarse block connections from tentative connections
   % listed in 'CG.cellFaces'.
   cf      = double(CG.cellFaces(is_act(CG.cellFaces(:,2)),:));
   cf(:,2) = renum(cf(:,2));

   % Determine number of half-faces connected to each coarse block.
   blk_siz = accumarray(cf(:,1), 1);

   % Build connections (at least topologically) for each coarse block.
   ix = @(n) [        repmat((1 : n).', [n, 1]), ...
              reshape(repmat( 1 : n   , [n, 1]), [], 1)];
   IJ = arrayfun(ix, blk_siz, 'UniformOutput', false);
   IJ = bsxfun(@plus, vertcat(IJ{:}), ...
               rldecode(cumsum([0; blk_siz(1:end-1)]), blk_siz.^2));

   B = sparse(IJ(:,1)       , IJ(:,2), 1);
   C = sparse(1 : size(cf,1), cf(:,1), 1); % Coarse-scale press gradient
   D = sparse(1 : size(cf,1), cf(:,2), 1); % Coarse-scale flux continuity

   % Build final synthetic matrix.
   A = [B  ,                        C    ,      D    ; ...
        C.', sparse(size(C,2), size(C,2) + size(D,2)); ...
        D.', sparse(size(D,2), size(C,2) + size(D,2))];
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function ix = bdry_supports_flow(G, CG, bc)
   if isempty(bc),
      ix = [];
   else
      [nsub, sub] = subFaces(G, CG);
      f2c         = sparse(sub, 1, ...
                           rldecode((1 : double(CG.faces.num)) .', nsub));

      flow_support                 = false([CG.faces.num, 1]);
      flow_support(f2c([bc.face])) = true;

      ix = find(flow_support);
   end
end
