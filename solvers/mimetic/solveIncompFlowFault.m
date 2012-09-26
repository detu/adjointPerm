function [xr, xw] = solveIncompFlowFault(xr, xw, g, s, fluid, varargin)
%Solve incompressible flow problem (fluxes/pressures) with fault and
%  shale multipliers
%
% SYNOPSIS:
%   [xr, xw] = solveIncompFlow(xr, xw, G, S, fluid)
%   [xr, xw] = solveIncompFlow(xr, xw, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (wells, sources, and boundary conditions).
%
% REQUIRED PARAMETERS:
%   xr, xw - Reservoir and well solution structures either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'solveIncompFlow' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, S   - Grid and (mimetic) linear system data structures as defined by
%            function 'computeMimeticIP'. G has to have defined shales
%            and fault multipliers from calculateFaultTrans.m to include
%            these effect
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
%   Solver       - Which solver mode function 'solveIncompFlow' should
%                  employ in assembling and solving the block system of
%                  linear equations.
%                  String.  Default value: Solver = 'hybrid'.
%
%                  Supported values are:
%                    - 'hybrid' --
%                         Assemble and solve hybrid system for interface
%                         pressures.  System is eventually solved by Schur
%                         complement reduction and back substitution.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','hybrid') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
%
%                    - 'mixed' --
%                         Assemble and solve a hybrid system for interface
%                         pressures, cell pressures and interface fluxes.
%                         System is eventually reduced to a mixed system as
%                         per function 'mixedSymm'.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','mixed') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
%
%                    - 'tpfa' --
%                         Assemble and solve a cell-centred system for cell
%                         pressures.  Interface fluxes recovered through
%                         back substitution.
%
%                         The system 'S' must in this case be assembled by
%                         passing option pair ('Type','mixed') or option
%                         pair ('Type','comp_hybrid') to function
%                         'computeMimeticIP'.
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
%                  the caller of function 'solveIncompFlow'.
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
%           'solveIncompFlow:DrivingForce:Missing'
%
%
% SEE ALSO:
%   computeMimeticIP, addBC, addSource, addWell, initSimpleFluid
%   initResSol, initWellSol, solveIncompFlowMS, calculateFaultTrans

%{
#COPYRIGHT#
%}

% $Id: solveIncompFlow.m 2342 2009-06-07 17:22:40Z bska $

   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'Solver',       'hybrid',         ...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false,            ...
                'condition_number', true);
   opt = merge_options(opt, varargin{:});

   g_vec   = gravity();
   no_grav = ~(norm(g_vec(1 : size(g.nodes.coords,2))) > 0);
   if all([isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.wells), no_grav]),
      warning(id('DrivingForce:Missing'),                      ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   else
      %add degre of freedom on fault faces, treat faults and shales at same
      %footing
      faultfaces = [];
      faultcellfaces = [];
      faulttrans = [];
      if(isfield(g,'shalefaces'))
          faultfaces = [faultfaces;g.shalefaces];
          faultcellfaces = [faultcellfaces;g.shalecellfaces];
          faulttrans = [faulttrans;g.shaletrans];
      end
      if(isfield(g,'faultfaces'))
          faultfaces = [faultfaces;g.faultfaces];
          faultcellfaces = [faultcellfaces;g.faultcellfaces];
          faulttrans = [faulttrans;g.faulttrans];
      end
      
      
      numf = ones(g.faces.num,1);
      numf(faultfaces) = repmat(2,numel(faultfaces),1);
      facenums = cumsum(numf);
      cellfaces =  facenums(g.cellFaces(:,1));
      if(~isempty(faultcellfaces))
        cellfaces(faultcellfaces(:,1))= cellfaces(faultcellfaces(:,2))-1;
      end
      
      
      %old cellfaces = g.cellFaces(:,1);
      numfaces=diff(g.cells.facePos);
      cellNo  = rldecode(1:g.cells.num, double(numfaces), 2) .';
      s.C     = sparse(1:numel(cellNo), cellNo, 1);
      s.D     = sparse(1:numel(cellNo), double(cellfaces), 1, ...
                        numel(cellfaces),facenums(end));
      %make the matrix elements correspongind to fault tramissibility
      mob   = fluid.Lt(xr);     % == total mobility for all cells
      omega = fluid.omega(xr);
      if(~isempty(faultcellfaces))
          trans = faulttrans;
          faultmob = 2./sum(1./mob(g.faces.neighbors(faultfaces,:)),2);
          trans = trans.*faultmob;
          
          trans = repmat(trans,4,1);
          trans(numel(trans)/2+1:end) = -trans(numel(trans)/2+1:end);
          trans = -1.0*trans;
          
          %ci = find(numf ==2);
          indl = cellfaces(faultcellfaces(:,1));
          indr = cellfaces(faultcellfaces(:,2));
          
          indll = [indl;indr;indl;indr];
          indrr = [indl;indr;indr;indl];
          
          s.E     = sparse(double(indll),double(indrr),trans,facenums(end),facenums(end));
      else
          s.E = sparse(facenums(end),facenums(end));
      end
      s.sizeB = repmat(size(g.cellFaces, 1), [1,2]);
      s.sizeC = size(s.C);
      s.sizeD = size(s.D);

      %[A, b, dF, dC] = build_system(xr, g, s, opt.wells, ...
      %                              opt.bc, opt.src, fluid, opt);
      % 
      %mob   = fluid.Lt(xr);     % == total mobility for all cells
       %omega = fluid.omega(xr);
      [Ar, br, dFr, dCr] = syscomp_res  (g, s, mob, omega, opt.bc, opt.src, opt);
      [Aw, bw, dFw, dCw] = syscomp_wells(g, opt.wells, mob, omega,          opt);

      A    = cell([1, 3]);           b    = cell([1, 3]);
      A{1} = blkdiag(Ar{1}, Aw{1});  b{1} = vertcat(br{1}, bw{1});
      A{2} = vertcat(Ar{2}, Aw{2});  b{2} =         br{2}        ;
      A{3} = blkdiag(Ar{3}, Aw{3});  b{3} = vertcat(br{3}, bw{2});

      %direchlet boundary conditions right side
      %b{3}(isdirech) = opt.bc.value(isdirech);
      
      
      
      % Finally, assemble a list of prescribed (contact) pressures, i.e., both
      % known face pressure and known bottom hole pressures.  These will be
      % eliminated from the linear system prior to calling one of the system
      % solvers, and entered directly into the solution component 'cp'.
      %
      dF   = [dFr; dFw];
      dC   = [dCr; dCw]; %%??? not used?                        
  
      numf = [numf;ones(size(dFw))];
      dF=rldecode(dF,numf);
      b{3}=rldecode(b{3},numf);
      %solver = @solve_hybrid;
      A{3} = A{3}(:,~dF);
      b{3} = b{3}(  ~dF);
      regul = ~any(dF);
      s.E = blkdiag(s.E,sparse(numel(dFw),numel(dFw)));
      %if(~isempty(faultcellfaces))
          E = s.E(~dF,~dF);
      %else
      %E=s.E;
      [x{1:5}] = schurComplementSymmFault(A{:}, E, b{:},...
                                     'Regularize', false, ...
                                     'LinSolve', opt.LinSolve,...
                                     'condition_number',opt.condition_number);
      lam = -ones(size(dF));                           
      lam(~dF) = x{3};
      x{3}     = lam;
      if(~isempty(faultcellfaces))
        fault_pressleft = x{3}(indl);
        fault_pressright = x{3}(indr);
      end
      if(length(x{3})>length(facenums))
        x{3} = [x{3}(facenums),x{3}(facenums(end)+1:end)];
      end
      if(~isempty(faultcellfaces))
        x{3}(faultfaces) = (fault_pressleft + fault_pressright);
      end
      
      
      %x        = solver(A, b);

      [xr, xw] = pack_solution(xr, xw, g, s, x{1:3}, opt);
      if(~isempty(faultcellfaces))
        xr.fault_pressleft  = fault_pressleft; 
        xr.fault_pressright  = fault_pressright;
      end
      if opt.MatrixOutput, 
          xr.A = x{4}; 
          xr.rhs=x{5};
      end
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function s = id(s)
   s = ['solveIncompFlow:', s];
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

function [xr, xw] = pack_solution(xr, xw, g, s, flux, pres, lam, opt)
   xr.cellPressure(:) = pres;
   %xr.facePressure(:) = lam (1 : s.sizeD(2)));
   xr.facePressure(:) = lam (1 : g.faces.num);
   xr.cellFlux(:)     = flux(1 : s.sizeB(1));
   xr.faceFlux(:)     = cellFlux2faceFlux(g, flux(1 : s.sizeB(1)));

   if ~isempty(opt.wells),
      nw  = numel(opt.wells);
      i_f = s.sizeB(1);
      i_p = g.faces.num;%s.sizeD(2);

      for k = 1 : nw,
         nperf = numel(opt.wells(k).cells);

         xw(k).flux     = -flux(i_f + 1 : i_f + nperf);
         xw(k).pressure =  lam (i_p + 1 : i_p +   1  );

         i_f = i_f + nperf;
         i_p = i_p +   1  ;
      end
   end
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_res(g, s, mob, omega, bc, src, opt)
   A = cell([1, 3]);

   mob = spdiags(s.C * mob, 0, s.sizeB(1), s.sizeB(2));

   if strcmpi(opt.Solver, 'hybrid'),
      if ~isfield(s, 'BI'),
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''hybrid'' is incompatible with ', ...
               'linear system type ''%s''.'], s.type);
      end
      A{1} = mob * s.BI;
   else
      if ~isfield(s, 'B'),
         error(id('SolverMode:Inconsistent'), ...
              ['Solver mode ''%s'' is incompatible with ', ...
               'linear system type ''%s''.'], opt.Solver, s.type);
      end
      A{1} = mob \ s.B;
   end

   A{2} = s.C;
   A{3} = s.D;


   [b{1:3}, grav, dF, dC] = computePressureRHS(g, omega, bc, src);
   b{1} = b{1} + grav;  % Add (mimetic) gravity contributions.
end

%--------------------------------------------------------------------------

function [A, b, dF, dC] = syscomp_wells(G, W, mob, omega, opt)
   % Assume empty well structure by default...
   A  = cell([1, 3]);
   b  = cell([1, 2]);
   dF = logical([]);
   dC = [];

   if ~isempty(W),
      % but fill in values when there nevertheless are some...
      nperf = cellfun(@numel, { W.cells });
      wc    = vertcat(W.cells);   n = numel(wc);   i = 1 : n;
      nW    = numel(W);

      % Diagonal transmissibility matrix (inner product).
      v = mob(wc) .* vertcat(W.WI);
      if ~strcmpi(opt.Solver, 'hybrid'), v = 1 ./ v; end
      A{1} = sparse(i, i, v, n, n); % == spdiags(v,0,n,n), but more direct.

      % Connection matrices {2} -> C and {3} -> D for all wells.
      A{2} = sparse(i, wc                      , 1, n, G.cells.num);
      A{3} = sparse(i, rldecode(1:nW, nperf, 2), 1, n, nW         );

      % Which (and what) are the prescribed well bottom-hole pressures?
      dF   = strcmpi('bhp', { W.type } .');
      dC   = reshape([ W(dF).val ], [], 1);

      % Form linsys rhs contributions, {1} -> pressure, {2} -> rate.
      b{1} = -vertcat(W.val);  b{1}(~dF) = 0;  % Remove rates
      b{2} = -vertcat(W.val);  b{2}( dF) = 0;  % Remove pressures

      % Expand well pressure rhs to each perforation, adjust for gravity
      % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
      dp   = norm(gravity()) * vertcat(W.dZ) .* omega(wc);
      b{1} = rldecode(b{1}, nperf) - dp;
   end
end

%--------------------------------------------------------------------------

function nF = neumann_faces(g, dF, opt)
% Determine the 'Neumann faces' (nF) (i.e., the faces for which
% face/contact pressures will be determined in the 'mixed' and 'TPFA'
% cases) of the model (g) according to the following rules:
%
%   - No internal faces are counted amongst the 'Neumann faces'.
%   - An external face is a 'Neumann face' unless there is a prescribed
%     pressure value (i.e., a Dirichlet condition) associated to the face.
%   - Wells controlled by rate-constraints are treated as Neumann contacts.
%
   nF                                 = false([g.faces.num, 1]); % No int.
   nF(any(g.faces.neighbors == 0, 2)) = true;  % All external faces ...
   nF(dF(1 : g.faces.num))            = false; % ... except Dirichlet cond.

   if ~isempty(opt.wells),
      % Additionally include all rate-constrained wells.
      nF = [nF; strcmpi('rate', { opt.wells.type } .')];
   end
end

%--------------------------------------------------------------------------

function Do = oriented_mapping(s, orient, opt)
% 'Do' maps face fluxes to half-face fluxes.  This matrix is used to form
% the reduced, mixed system of linear equtions in linear system solver
% functions 'mixedSymm' and 'tpfSymm'.
%
   if ~isempty(opt.wells),
      n  = sum(cellfun(@numel, { opt.wells.cells }));
      dw = speye(n);
      ow = ones([n, 1]);
   else
      dw = [];
      ow = [];
   end
   nf = numel(orient) + numel(ow);
   Do = spdiags([orient; ow], 0, nf, nf) * blkdiag(s.D, dw);
end
