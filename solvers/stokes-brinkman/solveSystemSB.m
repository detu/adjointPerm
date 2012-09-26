function [S,sol] = solveSystemSB(S,G, Dofs, varargin)
% solveSystemSB - Solves system S with boundary conditions BC
%
% SYNOPSIS:
%   [S,sol] = solveSystemSB(S,G, Dofs,Coord)
%   [S,sol] = solveSystemSB(S,G, Dofs,Coord, 'pn', pv, ...)
%
% DESCRIPTION:
%   Solves the 2-D or 3-D symmetric, mixed system of linear equations
%
%       [B1   0   C1] [v1] [F1]      [B1   0   0   C1] [v1] [F1]
%       [0    B2  C2] [v2]=[F2] or   [0    B2  0   C2] [v2]=[F2]
%       [C1'  C2' PO] [-p] [Q ]      [0    0   B3  C3] [v3]=[F3]
%                                    [C1'  C2' C3' PO] [-p] [Q ]
%
%   with respect to the x-, y- (and z-)velocities 'v1', 'v2' (and 'v3')
%   defined in all DOFS and the pressure 'p' in each vertice point. 
%
% PARAMETERS:
%   S          - The mixed Stokes-Brinkman system
%
%   G          - The grid structure
%
%   Dofs       - The degrees-of-freedom structure
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   regValue   - The value for p in regNode
%  
%   regNode    - The node to set the pressure level   
%
%   bc         - Boundary condition structure as defined by function 'addBCSB'.
%                This structure accounts for all external boundary conditions
%                to the reservoir flow.  May be empty (i.e., bc = struct([]))
%                which is interpreted as all external no-flow conditions.
%
%   src        - Explicit source contributions as defined by function
%                'addSource'.  May be empty (i.e., src = struct([])) which is
%                interpreted as a reservoir model without explicit sources.
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%
% RETURNS:
%   S           - The updated system structure with boundary conditions implemented
%
%   sol         - A structure containing the velocities 'v1', 'v2' (and
%                 'v3') and the pressure 'p'
  
  error(nargchk(4, 12, nargin, 'struct'));

  opt = struct('regValue', 0, 'regNode', 1, 'bc', [], 'src', [],'LinSolve', @mldivide);
  
  opt = merge_options(opt, varargin{:});
  if all([isempty(opt.bc), isempty(opt.src)]),
    warning(id('DrivingForce:Missing'),                      ...
            ['No external driving forces present in model--', ...
             'state remains unchanged.\n']);
  end
  
  dim = numel(G.cartDims);
 
  if isstruct(opt.src)
    S = setRHS_Q(S, Dofs.Pdofs, G, opt.src, 'src');
  end

  if ~isstruct(opt.bc)
    % Default: no-flow boundary conditions
    % NB: already applied to S in makeSystemSB for element = 'TH'
    if  ~strcmp(S.element,'TH')
       ind     = any(G.faces.neighbors==0,2);
       bcfaces = find(ind);
       BC      = addBCSB([],bcfaces, 'velocity_n', ...
                        repmat(0,numel(bcfaces),1), G, Dofs);
       S       = setBCSB_vn(S, BC);
    end
    BC.p.unique = [opt.regNode,opt.regValue];
    S           = setBCSB_p_unique(S,BC);
  else
     S = setBCSB(S,opt.bc, Dofs,G);
  end
  
  sb1 = size(S.B1); sb11 = sb1(1); sb12 = sb1(2);
  sb2 = size(S.B2); sb21 = sb2(1); sb22 = sb2(2);
  if dim==3
    sb3 = size(S.B3); sb31 = sb3(1); sb32 = sb3(2);
  end
  
  beta = 1./(max(max(abs([S.B1,S.B2])))/max(max(abs([S.C1,S.C2]))));

if isfield(S, 'D1')
  if dim==2
    B = [S.B1*beta         sparse(sb21,sb22) S.C1     ;...
         sparse(sb11,sb12) S.B2*beta         S.C2     ;...
         S.D1              S.D2              S.P0/beta];
    F = [S.F1    ;...
         S.F2    ;...
         S.Q/beta];
  elseif dim==3
    B = [S.B1*beta         sparse(sb21,sb22) sparse(sb31,sb32) S.C1   ;...
         sparse(sb11,sb12) S.B2*beta         sparse(sb31,sb32) S.C2   ;...
         sparse(sb11,sb12) sparse(sb21,sb22) S.B3*beta         S.C3   ;...
         S.D1              S.D2              S.D3              S.P0/beta];
    F = [S.F1    ;...
         S.F2    ;...
         S.F3    ;...
         S.Q/beta];
  end

else
  if dim==2
    B = [S.B1*beta         sparse(sb21,sb22) S.C1     ;...
         sparse(sb11,sb12) S.B2*beta         S.C2     ;...
         S.C1'             S.C2'             S.P0/beta];
    F = [S.F1    ;...
         S.F2    ;...
         S.Q/beta];
  elseif dim==3
    B = [S.B1*beta         sparse(sb21,sb22) sparse(sb31,sb32) S.C1   ;...
         sparse(sb11,sb12) S.B2*beta         sparse(sb31,sb32) S.C2   ;...
         sparse(sb11,sb12) sparse(sb21,sb22) S.B3*beta         S.C3   ;...
         S.C1'             S.C2'             S.C3'             S.P0/beta];
    F = [S.F1    ;...
         S.F2    ;...
         S.F3    ;...
         S.Q/beta];
  end
end

solution = opt.LinSolve(B, F);

  sol.v1 = solution(                1:  Dofs.numVdofs)*beta; 
  sol.v2 = solution(  Dofs.numVdofs+1:2*Dofs.numVdofs)*beta;
  if dim==2
    sol.p  = -solution(2*Dofs.numVdofs+1:end);
  elseif dim==3
    sol.v3 =  solution(2*Dofs.numVdofs+1:3*Dofs.numVdofs)*beta;
    sol.p  = -solution(3*Dofs.numVdofs+1:end);
  end
 