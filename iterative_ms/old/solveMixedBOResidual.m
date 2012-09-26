function [corrFunc, resSol, wellSol] = solveMixedBOResidual(cellInx, cellInxOverlap, ...
                                                             resSol, wellSol, ...
                                                             G, rock, S, fluid, ...
                                                             p0, dt, varargin)
% solveBOResidual  -- Version of solveBlackOilWellSystem for solving a
%                      local problem based on cellInx using mixed repr.
%
% SYNOPSIS:
%   [resSol, wellSol] = solveBOSubSystem(cellInx, resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt)
%   [resSol, wellSol] = solveBOSubSystem(cellInx, resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt,  ...
%                                               'pn1', pv1, ...)
%
% PARAMETERS:
%   cellInx - Indices to cells defining sub domain
%  
%   resSol  - Reservoir solution structure for full fine-grid from previous time step (or
%             previous iteration of successive substitution algorithm).
%
%   wellSol - Well solution structure for all wells from previous time step (or previous
%             iteration of successive substitution algorithm).
%
%   G       - Grid structure as described in SAMSIMGRID.
%
%   rock    - Rock structure.
%
%   S       - Linear system structure as defined by function
%             assembleMimeticSystem.
%
%   W       - Well linear system structure as defined by function
%             assembleWellSystem.
%
%   PVTTAB  - PVT tables.
%
%   p0      - Vector, length G.cells.num, of cell pressures at previous
%             time step.
%
%   dt      - Time step size.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - bc  -- Boundary condtion structure as defined by function
%                        'addBC'.  This structure accounts for all external
%                        boundary contributions to the reservoir flow.
%                        Default value: bc = [] meaning all external
%                        no-flow (homogeneous Neumann) conditions.
%
%               - src -- Explicit source contributions as defined by
%                        function 'addSource'.
%                        Default value: src = [] meaning no explicit
%                        sources exist in the model.
%
%               - LinSolve -- Solver for the resulting linear systems.
%                             Handle to a function supporting the syntax
%                                 x = LinSolve(A,b)
%                             for solving the system Ax=b of linear
%                             equations.
%                             Default value: LinSolve = @mldivide
%                             (backslash).
%
% RETURNS:
%   resSol  - Updated sparse reservoir solution structure.
%
%   wellSol - Updated 'sparse' well solution structure.


%--------------------------------------------------------------------------
%

opt = struct('bc', [], 'src', [], 'wells', [], ...
             'LinSolve', @mldivide, 'overlap', 0);
opt = merge_options(opt, varargin{:});
W = opt.wells;

%
subCellsOverlap = CG.cells.subCells;
if overlap > 0,
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
end
%% New approach
for k = 1:CG.cells.num
    [G1, S1, rock1, W1, bc1, indB, indD] = subProblem(G, S, rock, indCO, ...
                       'wells', W , 'bc', opt.bc, 'src', opt.src, ...
                       'subDomain', indC);
    [rs, ws] = subSol(resSol, )               

for k = 1:CG.cells.num
    indC  = CG.cells.subcells(:, k);
    indCO = subCellsOverlap(:, k);
    [resSol_k, wellSol_k, its] = solveRes(resSol, wellSol, G, rock, S, fluid, ...
                                     p0, dt, indC, indCO, varargin);
    corrFunc{k}.resSol  = resSol_k;
    corrFunc{k}.wellSol = wellSol_k;
    corrFunc{k}.its     = its;
    %XXX might have possibility to update resSol/wellSol inside loop 
end

if nargout > 1
    [resSol, wellSol] = addSols(resSol, wellSol, corrFunc);
end

%--------------------------------------------------------------------------
                                   
function [resSol_k, wellSol_k] = solveRes(resSol, wellSol, G, rock, S, fluid, ...
                                     p0, dt, indC, indCO, varargin)
resSol_k    = initSparseResSol(G);
wellSol_k   = initWellSol(W);
resSolTemp  = resSol;
wellSolTemp = wellSol; 

tol   = 1e-5,
pDiff = max(resSol.cellPressure) - min(resSol.cellPressure);
pTol  = tol;

vMax = max(abs(resSol.faceFlux));
vTol = tol;

pCrit = inf; vCrit = inf;

[Br, Cr, Dr, Pr, fr, gr, hr] = hybLinSys(resSolTemp, G, S, rock, fluid, p0, dt, ...
                                      'bc', opt.bc, 'src', opt.src);
[Bw, Cw, Dw, fw, gw, hw]     = hybLinWellSys(resSolTemp, wellSolTemp, G, S, ...
                                         W, rock, fluid);
                                     
%[Dred, Do, h] = mixedMappings(G, W, D, Dw, h, hw, opt.bc);
% Hybrid matrices
B = blkdiag(Br, Bw);
C = [Cr ; Cw];  
D = blkdiag(Dr, Dw);

indB  = logical( C * indC   );
indBO = logical( C * indCO  );
indD  = logical( D' * indB  );
indDO = logical( D' * indBO );

[locBnd, globPBnd, globNBnd] = getBoundaryType(G, W, opt.bc)

[]

subC = C(indBO, indCO);
subD = D(indBO, indDO);

[fluxCur, pCur, lamNCur]    = unpackSol(G, W, resSolTemp, wellSolTemp, opt.bc);

its = 0;
while ~done    
   subB = Do' * B(indBO, indBO) * Do;
   subP = Pr(indCO, indCO);
    
   f = [fr; fw] .* indB;  subf = Do' * f(indBO);
   g = gr .* indC;        subg = g(indCO);
   %h = [hr; hw] .* indD;  subh = h(indDO);
    
   fRes  = subf  - subB * subFluxCur + subC * subPCur - subD * sublamNCur;  
   gRes  = subg  - subC' * subFluxCur + subP * subPCur;
   hNRes = subhN - subD' * subFluxCur; 
   
   if (norm(fRes)/norm(pDiff)) < pTol ||
       done = true;
   else
       [flux, p, lamN] = solveMixed(B, C, D, fRes, gRes, hRes); 
       % Update fields
       subFluxCur  = subFluxCur + flux;
       subPCur     = subPCur + p;
       subLamNCur  = lamNCur + lamN;
       % Update system for next iteration
       fluxCur(indD)  = subFluxCur;
       pCur(indCO)    = subPCur;
       lamNInx()
       lamNCur(indDN) = subLamNCur(indHGlob);       
       
   [flux, p, lamN] = solveMixed
   fr = 
   [flux, p, lamN] = solveMixedLinSys(B, C, D, P, f, g, h, Do);
   
                                         
    % P does not change
    [B, C, D, f, g, h] = mixedSys(B, Bw, C, Cw, D, Dw, f, fw, g, gw, h, hw, opt.bc);
    
 
    
    [flux, p, lamN]    = unpackSol(G, W, resSolTemp, wellSolTemp, opt.bc);
    
    % Residuals:
    f = f - B*flux  + C*p - D*lamN;
    g = g - C'*flux - P*p;
    h = h - D'*flux;
    
    

% 


B = blkdiag(B, Bw);
C = [C ; Cw];
[D, Do, h] = mixedMappings(G, W, D, Dw, h, hw, opt.bc);

f = [f; fw]; 
g = g + gw;





%-----------------------------------------------------------------------
% Get updated pressure and saturation dependent quantities -------------
[c, rho, mu, u] = fluid.pvt(resSol.cellPressure, resSol.z);
s      = bsxfun(@rdivide, u, sum(u,2));
mob    = fluid.relperm(s) ./ mu;   
Lti    = 1./sum(mob, 2);   


% System matrices
nc     = G.cells.num;
cf     = G.cellFaces(:,1);
ncf    = numel(cf);
cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';

B  = spdiags(Lti(cellNo), ncf, ncf) * S.B;
C  = sparse((1:ncf)', cellNo, 1, ncf, nc);
D  = sparse((1:ncf)', double(cf), 1, ncf, G.faces.num);
accum  = sum(c, 2) .* poreVolume(G, rock) ./ dt;
P      = spdiags(-accum,0,nc,nc);

% System RHS
grav   = gravity(); grav = grav(:);
hg     = ( G.faces.centroids(cf, :) - G.cells.centroids(cellNo, :) ) * grav;
rhoTot = sum(rho.*mob ,2).*Lti;
f      = hg .* rhoTot(cellNo);

[f, g, h, gp, dF, dC] = computePressureRHS(G, rhoTot, opt.bc, opt.src);
a1   =  sum( c .* mob, 2) .* Lti;
a2   =  sum( c  .* mob .* bsxfun(@minus, 2*rhoTot, rho), 2) .* Lti;  % sjekk fortegn???
a3   =  rhoTot .* sum( c .* mob .* bsxfun(@minus, rhoTot, rho));     % sjekk fortegn???

g = g + g_comp + g_grav + S.RHS.volume_discrepancy;

f = [f + gp; fW];
h = [h     ; hW];


[f, g, h] = linHybRHS()
   % Compute gravity terms.
   grav   = gravity();
   grav   = reshape(grav(1 : size(G.nodes.coords,2)), 1, []);

   nc     = G.cells.num;
   cf     = G.cellFaces(:,1);
   cellNo = rldecode(1 : nc, double(G.cells.numFaces), 2) .';
   sgn    = 2*double(G.faces.neighbors(cf, 1) == cellNo) - 1;

   % cfn  = cell-face normal = *outward* face normals on each cell.
   cfn    = bsxfun(@times, G.faces.normals(cf,:), sgn);

   % ng  == n' *     g for all cellfaces.
   % nKg == n' * K * g for all cellfaces.
   [K, row, col] = permTensor(rock, size(G.nodes.coords,2));
   ng     = cfn * reshape(grav, [], 1);
   nKg    = sum(cfn(:,row) .* bsxfun(@times, K(cellNo,:), grav(col)), 2);

   %-----------------------------------------------------------------------
   % Get updated pressure and saturation dependent quantities -------------
   %
   [c, rho, mu, u] = fluid.pvt(resSol.cellPressure, resSol.z);
   s      = bsxfun(@rdivide, u, sum(u,2));
   mob    = bsxfun(@rdivide, fluid.relperm(s), mu);

   % Parameters used in pressure equation.
   g2     = grav * grav.';      % Recall: 'grav' is a *ROW* vector.
   Lt     = sum(       mob                              , 2)      ;
   ct     = sum( s  .*  c                               , 2)      ;
   omega  = sum(rho .* mob                              , 2) ./ Lt;
   zeta   = sum( c  .* mob                              , 2) ./ Lt;
   beta   = sum( c  .* mob .* bsxfun(@minus, omega, rho), 2)      ;

   %% Mimetic system with bc
   BIV    = bsxfun(@times, ng, beta(cellNo)) + ...
            zeta(cellNo) .* (resSol.cellFlux + omega(cellNo).*nKg);
   BIV    = sparse(1 : numel(cellNo), cellNo, BIV, numel(cellNo), nc);

   accum  = ct .* poreVolume(G, rock) ./ dt;
   P      = sparse(1 : nc, 1 : nc, -accum);  % == spdiags(-accum,0,nc,nc)
   g_comp = p0 .* accum;

   g_grav = - G.cells.volumes .* beta .* omega .* g2;

   %% Wells

   % Assume empty well structure by default...
   [BIW, BIVW, CW, DW, fW, hW] = deal([]);
   dFw = logical([]);
   dCw = [];

   if ~isempty(opt.wells),
      % but fill in values when there nevertheless are some...
      nperf = cellfun(@numel, { opt.wells.cells });
      wc    = vertcat(opt.wells.cells);   n = numel(wc);   i = 1 : n;
      nW    = numel(opt.wells);

      % Diagonal transmissibility matrix (inner product).
      BIW   = sparse(i, i, Lt(wc) .* vertcat(opt.wells.WI), n, n);

      % Connection matrices for all wells.
      CW    = sparse(i, wc                      , 1, n, nc);
      DW    = sparse(i, rldecode(1:nW, nperf, 2), 1, n, nW);

      % Which (and what) are the prescribed well bottom-hole pressures?
      dFw   = strcmpi('bhp', { opt.wells.type } .');
      dCw   = reshape([ opt.wells(dFw).val ], [], 1);

      % Form linsys rhs contributions, fW -> pressure, hW -> rate.
      fW    = -vertcat(opt.wells.val);  fW(~dFw) = 0;  % Remove rates
      hW    = -vertcat(opt.wells.val);  hW( dFw) = 0;  % Remove pressures

      % Expand well pressure rhs to each perforation, adjust for gravity
      % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
      dp    = norm(gravity()) * vertcat(opt.wells.dZ) .* omega(wc);
      fW    = rldecode(fW, nperf) - dp;

      BIVW  = sparse(i, wc, -zeta(wc) .* vertcat(wellSol.flux), n, nc);
      clear i n
   end

   %-----------------------------------------------------------------------
   % Build hybrid system components ---------------------------------------
   %
   n   = numel(cellNo);         i = 1 : n;
   C   = sparse(i, cellNo, 1);
   D   = sparse(i, double(G.cellFaces(:,1)), 1, n, G.faces.num);

   BI  = blkdiag(sparse(i, i, Lt(cellNo), n, n) * S.BI, BIW);     clear i n
   BIV = [       BIV ; BIVW];
   C   = [        C  ;  CW ];
   D   = blkdiag( D  ,  DW );

   [f, g, h, gp, dF, dC] = computePressureRHS(G, omega, opt.bc, opt.src);
   g = g + g_comp + g_grav + S.RHS.volume_discrepancy;

   f = [f + gp; fW];
   h = [h     ; hW];

   %--------------------------------------------------------------------------
   % residual RHS --------------------------------
   %
   v   = resSol.cellFlux;
   p   = resSol.cellPressure;
   sD  = sum(D)';
   Dp  = spdiags(1./sD, 0, numel(sD), numel(sD));

   pi = Dp*D'*(- BI\v + C*p);
   %pi  = D \ (- BI\v + C*p);
   %fr  = - BI\v + C*p - D*pi;
   gr  = g - C'*v + BIV'*(BI\v) + P*p; %(p-p0);
   hr  = h - D'*v;

   %--------------------------------------------------------------------------
   % Produce local indices --------------------------------
   %
   if ~islogical(cellInx)
      indC = logical( sparse(cellInx, 1, 1, G.cells.num, 1) );
   else
      indC = cellInx;
   end
   if ~islogical(cellInxOverlap)
      indCO = logical( sparse(cellInxOverlap, 1, 1, G.cells.num, 1) );
   else
      indCO = cellInxOverlap;
   end
   indB     = logical( C * indC );
   indD     = logical( D' * indB );

   indBO    = logical( C * indCO );
   indDO    = logical( D' * indBO );

   %--------------------------------------------------------------------------
   % Solve local hybrid system (unsymmetric) --------------------------------
   %
   %{
   if ~isempty(W),
      dFw = strcmp({ W.type } .', 'bhp');
      dCw = [ W(dFw).val ].';
   else
      dFw = logical([]);
      dCw = [];
   end
   %}
   dF = [dF; dFw];
   dC = [dC; dCw];

   % sub matrices:
   %indDD  = and(indD, dF);
   %----------------------------------
   indDnD = and(indDO, ~dF);
   %indDnD = indDO;

   fr = f - BI\v + C*p - D(:, ~dF)*pi(~dF);
   %fr = - BI\v + C*p - D*pi;
   %----------------------------------
   BIs  = BI(indBO, indBO);
   Cs   = C(indBO, indCO);
   %dFs  = dF(indD);
   Ds   = D(indBO, indDnD );
   BIVs = BIV(indBO, indCO);
   Ps   = P(indCO, indCO);
   fl   = fr.*indB; fs   = fl(indBO);
   gl   = gr.*indC; gs   = gl(indCO);
   %----------------------------------
   hs   = hr( indDnD );
   %hs   = zeros(size(hs));
   %----------------------------------

   flux    = sparse(size(S.B,2), 1);
   p       = sparse(G.cells.num, 1);
   lam     = sparse(size(D,2)  , 1);
   lam(dF) = dC;
   lam     = lam.*indD;

   [flux(indBO), p(indCO), lam(indDnD)] = schurComplement(BIs, Cs, Ds, BIVs, Ps, ...
                                                           fs, gs, hs,           ...
                                                          'LinSolve', opt.LinSolve);

   %--------------------------------------------------------------------------
   % Package solution in accessible form -------------------------------------
   %

   resSol_sub.cellPressure = p;
   resSol_sub.facePressure = lam (1 : G.faces.num);
   resSol_sub.cellFlux     = flux(1 : size(S.B,1));
   resSol_sub.faceFlux     = cellFlux2faceFlux(G, flux(1 : size(S.B,1)));

   if ~isempty(opt.wells),
      i_f = size(G.cellFaces, 1);
      i_p = G.faces.num;

      for k = 1 : nW,
         wellSol_sub(k).flux     = -flux(i_f + 1 : i_f + nperf(k)); %#ok
         wellSol_sub(k).pressure =  lam (i_p + 1 : i_p +    1    ); %#ok

         i_f = i_f + nperf(k);
         i_p = i_p +    1    ;
      end
   else
      wellSol_sub = struct([]);
   end


   
