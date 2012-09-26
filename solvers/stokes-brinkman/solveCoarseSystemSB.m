function [xr, CS]=solveCoarseSystemSB(G, CG, CS, Dofs, rock, fluid, part, element, varargin)
% Solve coarse (multiscale) system.
%
% SYNOPSIS:
%  [xr, CS] = solveCoarseSystemSB(G, CG, CS, Dofs, rock, fluid, part, element)
%  [xr, CS] = solveCoarseSystemSB(G, CG, CS, Dofs, rock, fluid, part, element,...
%                                 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (sources and boundary conditions).
%
% REQUIRED PARAMETERS:
%   G, CG   - Grid and coarse grid.
%            
%   CS      - Linear system structure on fine grid (S) and coarse grid
%            (CS) as defined by functions 'makeSystemSB' and
%            'generateCoarseSystemSB', respectively.
%
%   Dofs    - Degrees-of-freedom structure as defined by function 'findCartDofs'.
%
%   rock    - Rock objects containing the field 'perm'.
%
%   fluid   - Fluid object as defined by function 'initSingleFluid' also containing the 
%            field 'fluid.mu_eff'.
%
%   part    - Cell-to-block partition vector.
%
%   element - The name of a file containing element matrices. 
%             Default element = 'TH' gives file 'ElementMat2D'/'3D'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   bc     - Boundary condition structure as defined by function 'addBCsb'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = struct([]))
%            which is interpreted as all external no-flow conditions.
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
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - blockPressure-- Pressure values for all blocks in the coarse
%                            model 'CG'.
%          - blockFlux    -- Flux values for all blocks in the coarse
%                            model 'CG'.
%
% NOTE:
%   If there are no external influences, i.e., if both the structures
%   'bc', and 'src' are empty, then the input values 'xr' and 'CS' are
%   returned unchanged and a warning is printed in the command window.
%   This warning is printed with message ID 'solveCoarseSystem_sb:DrivingForce:Missing'
%   

   opt = struct('bc', [], 'src', [],'LinSolve', @mldivide);
   opt = merge_options(opt, varargin{:});

   if all([isempty(opt.bc), isempty(opt.src)]),
     warning(id('DrivingForce:Missing'),                      ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   else
     xr      = initResSol(G, 0.0); 
  [A, b, dF, dC, Psi] = build_system(G, CG, part, CS, Dofs, rock, fluid, element, ...
                                     opt);                                   
     CS.psi = Psi;
     x      = solve_hybrid(A, b, dF, opt.LinSolve);

     xr     = pack_solution(xr, G, CG, part, CS, Dofs, x{1:2});
   end
end
  
%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [A, b, dF, dC, Psi] = build_system(G, CG, p, CS, Dofs, rock, fluid,...
                                            element,opt)                                                    
         
   [nsub, sub] = get_subfaces(G, CG, CS);
   [I, J]      = coarsening_ops(G, CG, CS, sub, nsub);
   
   % Build system components B, C, D, f, g, and h for final hybrid system
   %
   %    [  B   C   D  ] [  v ]     [ f ]
   %    [  C'  0   0  ] [ -p ]  =  [ g ] .
   %    [  D'  0   0  ] [ cp ]     [ h ]
   %
   % The components are represented as one-dimensional cell arrays as
   %
   %    A = { B, C, D }
   %    b = { f, g, h }
   %
   % both of which are assembled separately from reservoir contributions
   % (syscomp_res).  
  
   %-----------------------------------------------------------------------
   % Assemble well sub matrices -------------------------------------------
   %
   [A, Psi] = syscomp_res(G, CG, CS, Dofs, p, opt);
   
   [b{1:3}, dF, dC] = sys_rhs(G, Dofs, CS, CG, Psi, I, J, opt.src, opt.bc);

   cellno = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2) .';
   
   % BDarcy is based on velocities, the coarse system is based on fluxes:
    actB    = CS.activeCellFaces; 
    actarea = cellFaceArea(CS, G, CG);
    area    = spdiags(1./actarea,0,size(actB,1),size(actB,1)); 

    % Construct BI = INV(A{1}) which is needed when solving hybrid system ---
   BI = sparse(numel(CS.activeCellFaces), numel(CS.activeCellFaces)); 
   for k = 1 : size(A{2},2), 
     cells = find(CG.cells.subCells(:,k)); 
     fhd   = unique(Dofs.HalfDofs(cells,:)); 
     if numel(G.cartDims)==2 
         fhd_psi = [fhd, fhd+Dofs.numHalfDofs];
     elseif numel(G.cartDims)==3   
         fhd_psi = [fhd, fhd+Dofs.numHalfDofs, fhd+2*Dofs.numHalfDofs];
     end
     range = find(A{2}(:,k)); 
     psi_k = Psi(fhd_psi, :);
     BDarcy = makeBDarcy(G, rock, fluid, Dofs, 'TH', cells);
     A{1}   = area* psi_k' * BDarcy * psi_k* area;
     BI(range, range) = inv( A{1}(range, range));
   end 
   A{1} = BI ; 
end 
 
%--------------------------------------------------------------------------

function [A, Psi] = syscomp_res(g, cg, cs, Dofs, p, opt)
% Form hybrid basis function matrices.
  
  Psi = hybrid_basis(g, cg, cs, Dofs);
  
  A    = cell([1, 3]);
  A{2} = cs.C(cs.activeCellFaces,        :      );
  A{3} = cs.D(cs.activeCellFaces, cs.activeFaces);

end
  
%--------------------------------------------------------------------------

function Psi = hybrid_basis(G, CG, CS, Dofs)
% Form hybrid basis function matrices Psi.
% Called from syscomp_res (only).
  
   % 1) Enumerate the (active) cell face unknowns according to the
   %    'cs.activeCellFaces' field.  As the basis functions are stored in
   %    the cs.activeFaces elements of the various 'cs.basis*' fields, we
   %    also create an alias for these coarse faces.
   %
   dim         = numel(G.cartDims);
   actB        = CS.activeCellFaces;
   actF        = CS.activeFaces;
   renum       = zeros([size(CG.cellFaces,1), 1]);
   renum(actB) = 1 : numel(actB);

   % 2) Build a block-to-face renumbering engine.  Specifically,
   %      f2hf_c(b,f) = j
   %    means that coarse face 'f' within coarse block 'b' corresponds to
   %    coarse flux unknown (coarse half-face) 'j' (when j > 0).
   %
   [nb, nf] = deal(double(CG.cells.num), double(CG.faces.num));
   f2hf_c   = sparse(double(CG.cellFaces(:,1)), ...
                     double(CG.cellFaces(:,2)), renum, nb, nf);

   % 3) Apply hybrid basis function splitting to each (active) flux basis
   %    function.  This process forms the SPARSE input triplet [i,j,v] from
   %    which the hybrid flux basis matrix Bv is finally formed.
   %
   %    Each element of the cell array { cs.basis{actF} } is a 'V' tuple as
   %    defined by function 'evalBasisFunc'.  We also need to renumber the
   %    block-to-face connections, so pass the static 'f2hf_c' matrix along
   %    with the dynamic 'V' tuples.
   %
   [i,j,v] = cellfun(@(x) split_for_th(x, f2hf_c, G, CG, Dofs), ...
                     { CS.basis{actF} }, 'UniformOutput', false);

   if dim==2
      Psi = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
      2*Dofs.numHalfDofs, numel(actB));
   elseif dim==3    
      Psi = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
      3*Dofs.numHalfDofs, numel(actB));
   end
end
   
%--------------------------------------------------------------------------

function [i, j, v] = split_for_th(c, hfno, G, CG, Dofs)
% Split a flux basis function tuple 'c' into a
% SPARSE matrix input triplet [i,j,v].  A hybrid basis function is split
% into one component for each of (at most) two coarse half-faces.  We
% recall that the input tuple, 'c', is a cell array of the form
%
%    1   2   3  4  5  6      for 2-D and
%   {i, v1, v2, f, b, n}
%
%    1   2   3   4  5  6  7   for 3-D
%   {i, v1, v2, v3, f, b, n}
%
% with symbols as defined in 'evalBasisFuncTH'.  Item 1 directly provides the
% desired output vector 'i' (no further processing required).
%
% Our task, then, is to convert item 2-4 into the output vector 'v', and
% to build output vector 'j' from scratch. The output 'v' is the
% simplest to define.  We simply negate the v-values (item 2-4)
% corresponding to the second block (b(2)) if applicable.  If
% numel(b)==1, then v==c{2} is returned unchanged.  The 'j' vector
% contains at most two distinct values, each corresponding to a specific
% coarse half-face identified by a pair (b,f) of coarse blocks and face.
  
  dim=numel(G.cartDims);
  
  % 1) Extract 'i' values (the fine half-dofs) directly from input tuple.
  if dim==2
      i = [c{1}; c{1}+Dofs.numHalfDofs];
  elseif dim==3
      i = [c{1}; c{1}+Dofs.numHalfDofs; c{1}+2*Dofs.numHalfDofs];
  end
  
  % 2) Build 'j' values from scratch by 
  %i - Defining pairs of coarse blocks and face.
  if dim==2
     b = reshape(double(c{5}), [], 1);
     f = repmat (double(c{4}), [numel(b), 1]);
  elseif dim==3
     b = reshape(double(c{6}), [], 1);
     f = repmat (double(c{5}), [numel(b), 1]);
  end
  
  %    ii  - Extracting coarse half-faces identified by these pairs.
  j = reshape(full(hfno(sub2ind(size(hfno), b, f))), [], 1);
  
  %    iii - Expanding these half-face indices a number of times equal to
  %          the number of entities (item 5) in each of the coarse blocks
  %          identified by b.
  if dim==2
      j = rldecode(j, reshape(c{6}, [], 1));
      j = [j; j];
  elseif dim==3
      j = rldecode(j, reshape(c{7}, [], 1));
      j = [j; j; j];
  end

  % 3) Negate b(2) values if numel(b) > 1.  Else somewhat expensive NOP.
  s = [1; -1];
  if (numel(b)==1)
    [nsub, sub] = subFaces(G, CG);
    sub_ix      = cumsum([0; nsub]);
    fface       = sub(sub_ix(f) + 1 : sub_ix(f + 1));
    cellf       = G.faces.neighbors(fface(1),:);

    iC = find(CG.cells.subCells(:,b));
    if(find(iC==cellf(1))) 
      cellf = cellf(1);
    else 
      cellf = cellf(2);
    end
    cellNo = rldecode(1:G.cells.num, double(G.cells.numFaces), 2).';
    hf     = find(cellNo==cellf&G.cellFaces(:,1)==fface(1));
    ftag   = G.cellFaces(hf,2);
    if (ftag==1 || ftag==3 || ftag==5)
      s = -1;
    else
      s = 1;
    end
  end
  
  if dim==2
     v1 = rldecode(s(1 : numel(b)), reshape(c{6}, [], 1)) .* c{2};
     v2 = rldecode(s(1 : numel(b)), reshape(c{6}, [], 1)) .* c{3};
     v  = [v1; v2];
  elseif dim==3
     v1 = rldecode(s(1 : numel(b)), reshape(c{7}, [], 1)) .* c{2};
     v2 = rldecode(s(1 : numel(b)), reshape(c{7}, [], 1)) .* c{3};
     v3 = rldecode(s(1 : numel(b)), reshape(c{7}, [], 1)) .* c{4};
     v  = [v1; v2; v3];
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

function [I, J] = coarsening_ops(G, CG, CS, sub, nsub)
% coarsening_ops -- Compute coarsening operators from fine to coarse model.
%
% SYNOPSIS:
%   [I, J] = coarsening_ops(G, CG, CS, sub, nsub)
%
% PARAMETERS:
%   G, CG     - Fine and coarse grid discretisation of the reservoir model,
%               respectively.
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

   I = CG.cells.subCells;
   J = sparse(sub, rldecode(CS.activeFaces(:), nsub), ...
              1, G.faces.num, CG.faces.num);
   J = J(:, CS.activeFaces);
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = sys_rhs(G, Dofs, CS, CG, Psi, I, J, src, bc)
% sys_rhs -- Evaluate coarse system right hand side contributions.
%
% SYNOPSIS:
%   [f, G, h, dF, dC] = sys_rhs(G, omega, Psi, I, J, src, bc)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   Psi   - Flux basis functions from reservoir blocks.
%
%   I, J  - Coarsening operators as defined by function 'coarsening_ops'.
%
%   bc    - Boundary condition structure as defined by function 'add3DBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Coarse-scale direct (block) right hand side contributions as
%             expected by the linear system solvers such as
%             'schurComplement'.
%
%   dF, dC  - Coarse-scale packed pressure condition structure.
%             The faces in 'dF' and the values in 'dC'.  May be used to
%             eliminate known face pressures from the linear system before
%             calling a system solver (e.g., 'schurComplement').
%
% SEE ALSO:
%   computePressureRHS.

   % Evaluate fine-scale system right hand side contributions.
   %
   [ff, gg, hh, dF, dC] = computePressureRHSSB(G, Dofs, bc, src);

   % Accumulate fine-scale boundary conditions to coarse scale.
   %
   %ff = Psi.' * ff; 
   gg = I .' * gg;
   hh = J .' * hh;

   % Account for coarse scale Dirichlet boundary conditions.
   %
   dC2 = zeros(size(dF));  dC2(dF) = dC;
   dF  = J.' * dF;
   dC  = J.' * dC2 ./ dF;
   dF  = logical(dF);
   dC  = dC(dF);
   
   %ff_a  = zeros(size(CG.cellFaces,1),1);
   ff_af = zeros(CG.faces.num,1);
   ff_af(CS.activeFaces(find(dF))) = -dC;
   ff_a  = ff_af(CG.cellFaces(:,2));
   ff    = ff_a(CS.activeCellFaces);
   
end

%--------------------------------------------------------------------------

function x = solve_hybrid(A, b, dF, LinSolve)

  regul = ~any(dF);
  A{3} = A{3}(:,~dF);
  b{3} = b{3}(  ~dF);
  
  [x{1:4}] = schurComplementSymm(A{:}, b{:},          ...
                                 'regularize', regul, 'LinSolve', LinSolve);
  lam(~dF) = x{3};
  x{3}     = lam;
end
   
function xr = pack_solution(xr, G, CG, p, CS, Dofs, flux, pres)
  
  [nsub, sub] = get_subfaces(G, CG, CS);
  actB = CS.activeCellFaces;
  actF = CS.activeFaces;
  dim  = numel(G.cartDims);
  xr   = rmfield(xr,'facePressure');
  xr   = rmfield(xr,'cellFlux');
  xr   = rmfield(xr,'s');
  %-----------------------------------------------------------------------
   %% Package reservoir (coarse) solution structure -----------------------
   %
  xr.blockPressure   = pres;
  xr.blockFlux       = flux(1 : numel(actB));
  xr.blockVel        = xr.blockFlux./cellFaceArea(CS, G, CG);
  xr.cellPressure(:) = pres(p);
  
  basis         = [];
  basis.basisx  = CS.psi(1:Dofs.numHalfDofs,:);
  basis.basisy  = CS.psi(Dofs.numHalfDofs+1:2*Dofs.numHalfDofs,:);
  if dim==3
    basis.basisz  = CS.psi(2*Dofs.numHalfDofs+1:3*Dofs.numHalfDofs,:);
  end
  xr            = convertMStoFSvel(xr, G, CG, Dofs, basis, CS.activeCellFaces);
  
end

%--------------------------------------------------------------------------

function Phi = res_basisP(G, CG, CS, p, mob)
% res_basisP -- Assemble pressure basis functions for reservoir blocks.
%
% SYNOPSIS:
%   Phi = res_basisP(G, CG, CS, p, mob)
%
% PARAMETERS:
%   G, CG - Fine-scale and coarse-scale discretisations (grids) of
%           reservoir model.
%
%   CS    - Coarse linear system structure as defined by function
%           'generateCoarseSystem'.
%
%   p     - Cell-to-coarse block partition mapping vector.
%
%   mob   - Vector of (current) total mobility values.  One scalar value
%           for each cell in the underlying fine-scale model.
%
% RETURNS:
%   Phi - Reservoir pressure basis functions represented as a sparse
%         matrix of size G.cells.num-by-nchc.

   lam  = accumarray(p, mob .* G.cells.volumes) ./ ...
          accumarray(p,        G.cells.volumes);
   dlam = CS.mobility(CS.activeCellFaces)       ./ ...
          lam(CG.cellFaces(CS.activeCellFaces, 1));

   nc  = numel(CS.activeCellFaces);
   Phi = CS.basisP(:, CS.activeCellFaces) * spdiags(dlam, 0, nc, nc);
end

%--------------------------------------------------------------------------

function [ff, gg, hh, dF, dC] = computePressureRHSSB(g, Dofs, bc, src)
% computePressureRHSSB -- Construct mimetic linear system RHS block
%                          contributions.
%
% SYNOPSIS:
%   [f, g, h, dF, dC] = computePressureRHSSB(G, bc, src)
%
% PARAMETERS:
%   G     - Grid data structure.
%
%   Dofs  - Degrees-of-freedom data structure.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions to
%           the reservoir flow.  May be empty (i.e., bc = struct([])) which
%           is interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = struct([])) which is
%           interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   f, g, h - Direct (block) right hand side contributions as expected by
%             the linear system solvers such as 'schurComplement'.
%
%   dF, dC  - Packed Dirichlet/pressure condition structure.  The faces in
%             'dF' and the values in 'dC'.  May be used to eliminate known
%             face pressures from the linear system before calling a system
%             solver (e.g., 'schurComplement').
%
% SEE ALSO:
%   addBC, addSource, pvt, assembleMimeticSystem, schurComplement.

% $Id: solveCoarseSystemSB.m 2924 2009-10-01 14:13:14Z bska $

   dim=numel(g.cartDims);
   ff = zeros(dim*Dofs.numHalfDofs,1,1);
   gg = zeros([g.cells.num, 1]);
   hh = zeros([g.faces.num, 1]);

   if ~isempty(src),
      assert (max([src.cell]) <= g.cells.num);

      ss = accumarray(src.cell, src.rate);
      ii = accumarray(src.cell, 1);
      gg(ii > 0) = gg(ii > 0) + ss(ii > 0);
   end

   dF = false([g.faces.num, 1]);
   dC = [];
   
   if ~isempty(bc),
      assert (max([bc.face]) <= g.faces.num);
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      isDir = bc.type == 2; %strcmp('pressure', bc.type);
      
      [face, ix] = sort(bc.face(isDir));
      dF(face)   = true;

      % Add Dirichlet/pressure conditions to RHS.
      dC = bc.value(isDir);
      dC = dC(ix);    % Order according to increasing face index.
   end

   assert (~any(dC < 0));  % Pressure conditions should always be non-neg.
end
