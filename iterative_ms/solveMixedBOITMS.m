function [xrMs, xwMs, itNo] = solveMixedBOITMS(xrMs, xwMs, G, CG, p, rock, S, CS, ...
                                         fluid, DT, varargin)
                                     
opt = struct('bc', [], 'src', [], 'wells', [], ...
             'LinSolve', @mldivide, 'overlap', 0, ...
             'tol', 1e-3, 'Verbose', true, ...
             'useCorrection', true, ...
             'singleCorrection', false, ...
             'maxIt', 50);
opt = merge_options(opt, varargin{:});                  

if opt.Verbose
    fprintf('******* Solving Pressure Equation by iterative MS: ******* \n');
end

[I, J, Phi, Psi] = restrictionOps(G, S, opt.wells, CG, p, CS);
LS = initMixedBO(G, opt.bc, opt.src, opt.wells, true);

p0   = xrMs.cellPressure;
[xrMs, xwMs, incn, resn, its] = ...
   solveMixedBOMS(xrMs, xwMs, G, CG, p, rock, S, CS, fluid, p0, DT, LS, ...
                  I, J, Phi, Psi, ...
                  'tol', opt.tol,'maxIt', opt.maxIt);
               
itNo      = zeros(1,CG.cells.num+1);
itNo(1)   = itNo(1)+its;
inc       = 0; %max(incn.flux,incn.pressure);
resp      = inf;
res       = max(resn.f,resn.g);
tolRes    = 0.5*opt.tol;
tolBOMS   = 0.5*opt.tol;
if opt.useCorrection
   it   = 0;
   while max(inc, res) > opt.tol && it < opt.maxIt % && res < 1.1*resp
      it = it+1;
        
      [corrFunc,its] = ...
         solveMixedBOResiduals(xrMs, xwMs, G, rock, S, CG, fluid, p0, DT, ...
                               LS, opt.overlap, 'tol', tolRes, ...
                               'maxIt', opt.maxIt, 'Verbose', false);
      itNo(2:end) = itNo(2:end)+its;

      LS = setupMixedBO(xrMs, xwMs, G, rock, S, fluid, p0, DT, LS);
      resn.f = norm(LS.f)/(max(xrMs.cellPressure) - min(xrMs.cellPressure));
      resn.g = norm(LS.g)/norm(xrMs.faceFlux);
%      clf; plotCellData(G,LS.C'*abs(LS.f),'edgecolor','k'); colorbar; view(2)
%      fprintf('   %0*d: Relative Norms: Res = %10.4e\n\n', 3, it, max(resn.f,resn.g));
      
      ax = caxis;
      for k = 1:numel(corrFunc)
         xrMs.cellPressure = xrMs.cellPressure + corrFunc(k).resSol.cellPressure;
         xrMs.faceFlux = xrMs.faceFlux + corrFunc(k).resSol.faceFlux;
         xrMs.facePressure = xrMs.facePressure + corrFunc(k).resSol.facePressure;
         xrMs.cellFlux = faceFlux2cellFlux(G, xrMs.faceFlux);
         LS = setupMixedBO(xrMs, xwMs, G, rock, S, fluid, p0, DT, LS);
         resn.f = norm(LS.f)/(max(xrMs.cellPressure) - min(xrMs.cellPressure));
         resn.g = norm(LS.g)/norm(xrMs.faceFlux);
%         clf; plotCellData(G,LS.C'*abs(LS.f),'edgecolor','k'); colorbar; view(2), caxis(ax);
%         fprintf('   %0*d: Relative Norms: Res = %10.4e\n\n', 3, it, max(resn.f,resn.g));
      end
      xrMs.cellFlux = faceFlux2cellFlux(G, xrMs.faceFlux);
            
      if opt.singleCorrection, return; end
      
      [xrMs, xwMs, incn, resn, its] = ...
         solveMixedBOMS(xrMs, xwMs, G, CG, p, rock, S, CS, fluid, p0, DT, ...
                        LS, I, J, Phi, Psi, ...
                        'tol', tolBOMS, 'maxIt', opt.maxIt);
      itNo(1) = itNo(1)+its;
      inc  = max(incn.flux,incn.pressure);
      resp = res;
      res  = max(resn.f,resn.g);
      if opt.Verbose
         fprintf('Done outer loop:\n')
         fprintf('   %0*d: Relative Norms: Inc = %10.4e  Res = %10.4e\n\n', 3, it, inc, res);
      end
      if inc < opt.tol, tolRes = max(tolRes/5, 0.1*opt.tol); end
   end
end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------
function [Psi, Phi] = resBasis(G, S, CG, CS)
    [BPsi, Phi] = basisMatrixMixed(G, CG, CS);
    Psi = S.BI * BPsi;
end


function [Psi, Phi] = wellBasis(G, W) 
   if isempty(W),
      Psi = sparse(size(G.cellFaces,1), 0);
      Phi = sparse(G.cells.num        , 0);
   else
       CS  = [ W.CS   ];
       [i,j,v] = cellfun(@(x) x{1:3}, vertcat(CS.basis), ...
                         'UniformOutput', false);
       Psi     = sparse(      vertcat(i{:}) , ...
                       cumsum(vertcat(j{:})), ...
                              vertcat(v{:}) , ...
                       size(G.cellFaces,1), numel(j));
      
       [i,j,v] = cellfun(@(x) x{1:3}, vertcat(CS.basisP), ...
                        'UniformOutput', false);
       Phi     = sparse(       vertcat(i{:}) , ...
                        cumsum(vertcat(j{:})), ...
                               vertcat(v{:}) , G.cells.num, numel(j));  
   end
end


function [nsub, sub] = get_subfaces(G, CG, CS)
% get_subfaces -- Extract fine-scale faces constituting all active coarse
%                 faces.
%
% SYNOPSIS:
%   [nsub, sub] = get_subfaces(G, CG, CS)
%
% PARAMETERS:
%   G, CG - Fine and coarse grid discretisation of the reservoir model,
%           respectively.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
% RETURNS:
%   nsub, sub - The return values from function 'subFaces' compacted to
%               only account for active coarse faces (i.e. coarse faces
%               across which there is a non-zero flux basis function and,
%               therefore a possibility of non-zero flux).

   [nsub, sub] = subFaces(G, CG);
   sub_ix      = cumsum([0; nsub]);
   sub         = sub(mcolon(sub_ix(CS.activeFaces) + 1, ...
                            sub_ix(CS.activeFaces + 1)));
   nsub        = nsub(CS.activeFaces);
end

%--------------------------------------------------------------------------

function [I, J] = coarsening_ops(G, CG, p, CS, sub, nsub)
% coarsening_ops -- Compute coarsening operators from fine to coarse model.
%
% SYNOPSIS:
%   [I, J] = coarsening_ops(G, CG, p, CS, sub, nsub)
%
% PARAMETERS:
%   G, CG, p  - Fine and coarse grid discretisation of the reservoir model,
%               as well as cell-to-coarse block partition vector.
%
%   CS        - Coarse linear system structure as defined by function
%               'generateCoarseSystem'.
%
%   sub, nsub - Packed representation of fine-scale faces constituting
%               (active) coarse faces.  Assumed to follow the conventions
%               of the return values from function 'subFaces'.
%
% RETURNS:
%  I - Coarse-to-fine cell scatter operator (such that I.' is the
%      fine-to-coarse block gather (sum) operator).
%
%  J - Coarse-to-fine face scatter operator for active coarse faces such
%      that J.' is the fine-to-coarse face gather (sum) operator.

   I = sparse(1 : G.cells.num, p, 1, G.cells.num, CG.cells.num);
   J = sparse(sub, rldecode(CS.activeFaces(:), nsub), ...
              1, G.faces.num, CG.faces.num);
   J = J(:, CS.activeFaces);
end

function [I, J, Phi, Psi] = restrictionOps(G, S, W, CG, part, CS)
[nsub, sub] = get_subfaces(G, CG, CS);
[I, J]      = coarsening_ops(G, CG, part, CS, sub, nsub);

[Psi, Phi]     = resBasis(G, S, CG, CS);
[Psi_w, Phi_w] =  wellBasis(G, W);

Psi = [Psi Psi_w];
Phi = [Phi Phi_w];
end