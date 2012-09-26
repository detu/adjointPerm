function [xr, xw] = incompTPFA(xr, xw, G, T, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%
% SYNOPSIS:
%   [xr, xw] = incompTPFA(xr, xw, G, T, fluid)
%   [xr, xw] = incompTPFA(xr, xw, G, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%   xr, xw - Reservoir and well solution structures either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompTPFA' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput - Whether or not to return the final system matrix 'A' to
%                  the caller of function 'incompTPFA'.
%                  Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%          - cellFlux     -- Outgoing flux across each local interface for
%                            all cells in the model.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'incompTPFA:DrivingForce:Missing'
%
% EXAMPLE:
%    G   = computeGeometry(cartGrid([3,3,5]));
%    f   = initSingleFluid();
%    rock.perm = rand(G.cells.num, 1)*darcy()/100;
%    bc  = pside([], G, 'LEFT', 1:G.cartDims(2), 1:G.cartDims(3),2);
%    src = addSource([], 1, 1);
%    W   = verticalWell([], G, rock, 1, G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'rate', 'Val', 1/day(), ...
%                       'InnerProduct', 'ip_tpf');
%    W   = verticalWell(W, G, rock, G.cartDims(1),   G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'bhp', ...
%                       'Val',  1*barsa(), 'InnerProduct', 'ip_tpf');
%    T   = computeTrans(G, rock);
%    xr  = initResSol (G, 10);
%    xw  = initWellSol(G, 10);
%    [xr,xw] = incompTPFA(xr, xw, G, T, f, 'bc',bc,'src',src,'wells',W,...
%                         'MatrixOutput',true);
%
%    plotCellData(G, xr.cellPressure);
%
% SEE ALSO:
%   computeTrans, addBC, addSource, addWell, initSingleFluid, initResSol,
%   initWellSol.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-26 15:02:21 +0100 (ma, 26 okt 2009) $
% $Revision: 3078 $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      false,...
                'condition_number',false);
            
   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   no_grav = ~(norm(g_vec(1 : size(G.nodes.coords,2))) > 0);
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), no_grav]),
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Setting up linear system...\t\t\t');
   ticif (opt.Verbose);

   % Preliminaries
   cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
   cf     = G.cellFaces(:,1);
   nf     = G.faces.num;
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   % Face transmissibility = harmonic average of half-transmissibilities
   totmob = fluid.Lt(xr);
   T      = T.*totmob(cellNo);
   ft     = 1./accumarray(cf, 1./T, [nf, 1]);

   % Identify internal faces
   i  = all(G.faces.neighbors ~= 0, 2);

   % Boundary conditions and source terms.
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(G, fluid.omega(xr), ...
                                                   opt.bc, opt.src);
   sgn = 2*(G.faces.neighbors(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]);

   rhs = accumarray(cellNo, -ft(cf).*(sgn.*fg(cf)+ff), [n, 1]) + ...
         [gg; zeros(nw, 1)]                                    + ...
         accumarray(cellNo, -hh(cf), [n, 1]);


   d  = zeros(G.cells.num, 1);

   % Wells --------------------------
   C    = cell (nw, 1);
   D    = zeros(nw, 1);
   W    = opt.wells;

   for k = 1 : nw,
      wc       = W(k).cells;
      nwc      = numel(wc);
      w        = k + nc;

      wi       = W(k).WI .* totmob(wc);

      dp       = norm(gravity()) * W(k).dZ*sum(fluid.rho .* W(k).compi, 2);
      d   (wc) = d   (wc) + wi;

      if     strcmpi(W(k).type, 'bhp'),
         ww=max(wi);
         %ww=1.0;
         rhs (w)  = rhs (w)  + ww*W(k).val; 
         rhs (wc) = rhs (wc) + wi.*(W(k).val + dp);
         C{k}     = -sparse(1, nc);
         D(k)     = ww;

      elseif strcmpi(W(k).type, 'rate'),
         rhs (w)  = rhs (w)  + W(k).val;
         rhs (wc) = rhs (wc) + wi.*dp;

         C{k}     =-sparse(ones(nwc, 1), wc, wi, 1, nc);
         D(k)     = sum(wi);

         rhs (w)  = rhs (w) - wi.'*dp;

      else
         error('Unsupported well type.');
      end
   end

   C = vertcat(C{:});
   D = spdiags(D, 0, nw, nw);
   %-----------------------------------------


   % Add up internal face transmissibilities plus Dirichlet pressure
   % faces for each cell.
   d  = d + ...
        accumarray(cellNo(dF(cf)), T(dF(cf)), [nc, 1]) +...
        accumarray(reshape(G.faces.neighbors(i,:), [], 1), ...
        repmat(ft(i), [2,1]),  [nc, 1]);

   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries.  Also, wells introduce
   % additional equations and unknowns.
   I  = [G.faces.neighbors(i,1); G.faces.neighbors(i,2); (1:nc)'];
   J  = [G.faces.neighbors(i,2); G.faces.neighbors(i,1); (1:nc)'];
   V  = [-ft(i); -ft(i); d];
   A  = sparse(double(I), double(J), V, nc, nc);
   A = [A, C'; C D];

   tocif(opt.Verbose);


   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Solving linear system...\t\t\t');
   ticif (opt.Verbose);

   if ~any(dF) && (isempty(W) || ~any(strcmpi({ W.type }, 'bhp'))),
       if(A(1)>0)
        A(1) = 2*A(1);
       else
           error('A(1) not > 0')
       end
           
   end
   if opt.condition_number,
      %tic;
      disp('***************************************');
      disp(['Conditon number is :   ', num2str(condest(A))]);
      disp('***************************************');
      %toc;
   end
   p = opt.LinSolve(A, rhs);

   tocif(opt.Verbose);

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t');
   ticif (opt.Verbose);

   % Reconstruct face pressures and fluxes.
   fpress     =  ...
          accumarray(G.cellFaces(:,1), (p(cellNo)+grav).*T, [G.faces.num,1])./ ...
          accumarray(G.cellFaces(:,1), T, [G.faces.num,1]);


   % Neumann faces
   b         = any(G.faces.neighbors==0, 2);
   fpress(b) = fpress(b) - hh(b)./ft(b);


   % Contribution from gravity
   %fg         = accumarray(cf, grav.*sgn, [nf, 1]);
   %fpress(~i) = fpress(~i) + fg(~i);

   % Dirichlet faces
   fpress(dF) = dC;


   % Sign for boundary faces
   sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   ni   = G.faces.neighbors(i,:);
   flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nf, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nf, 1]);
   flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );
   %flux = -sgn.*ft((fpress(~i)-p(c)-grav));
   xr.cellPressure(:) = p(1:nc);
   xr.faceFlux(:)     = flux;
   xr.facePressure    = fpress;
   xr.cellFlux        = faceFlux2cellFlux(G, flux);

   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(fluid.rho.*W(k).compi, 2);
      xw(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      xw(k).pressure = p(nc+k);
   end

   if opt.MatrixOutput, 
       xr.A = A; 
       xr.rhs=rhs;
   end

   tocif(opt.Verbose);
end
