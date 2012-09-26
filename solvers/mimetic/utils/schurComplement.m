function [flux, press, lam] = schurComplement(BI, C, D, BIV, P, f, g, h, ...
                                              varargin)
%Solve non-symmetric system of linear eqns using Schur Complement analysis.
%
% SYNOPSIS:
%   [flux, press, lam] = schurComplement(BI, C, D, BIV, P, f, g, h)
%   [flux, press, lam] = schurComplement(BI, C, D, BIV, P, f, g, h, ...
%                                        'pn1', pv1, ...)
%
% DESCRIPTION:
%   Solves the unsymmetric, hybrid (block) system of linear equations
%
%       [    B    C   D  ] [  v ]     [ f ]
%       [  C'-V'  P   0  ] [ -p ]  =  [ g ]
%       [    D'   0   0  ] [ cp ]     [ h ]
%
%   with respect to the half-face (half-contact) flux 'v', the cell
%   pressure 'p', and the face (contact) pressure 'cp'.  The system is
%   solved using a Schur complement reduction to an unsymmetric system
%   of linear equations
%
%        S cp = R    (*)
%
%   from which the contact pressure 'cp' is recovered.  Then, using back
%   substitution, the cell pressure and half-contact fluxes are
%   recovered.
%
%   Note, however, that the function assumes that any faces having
%   prescribed pressure values (such as Dirichlet boundary condtions or
%   wells constrained by BHP targets) have already been eliminated from
%   the blocks 'D' and 'h'.
%
% PARAMETERS:
%   BI      - Inverse of the matrix B in the block system.  Assumed SPD.
%   C       - Block 'C' of the block system.
%   D       - Block 'D' of the block system.
%   BIV     - Matrix product BI*V, with 'V' of the block system.
%   P       - Block 'P' of the block system.  Assumed diagonal.
%   f       - Block 'f' of the block system right hand side.
%   g       - Block 'g' of the block system right hand side.
%   h       - Block 'h' of the block system right hand side.
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               LinSolve -- Function handle to a solver for the system (*)
%                           above.  Assumed to support the syntax
%
%                                    x = LinSolve(A, b)
%
%                           to solve a system Ax=b of linear equations.
%                           Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   flux  - The half-contact fluxes 'v'.
%   press - The cell pressure 'p'.
%   lam   - The contact pressure 'cp'.

%{
#COPYRIGHT#
%}

% $Id: schurComplement.m 1953 2009-03-31 10:54:12Z bska $

opt = struct('LinSolve', @mldivide);
opt = merge_options(opt, varargin{:});

%--------------------------------------------------------------------------
% Schur complement analysis to contact pressure system --------------------
%
ncell = numel(g);               % Number of cells in grid

BIDf   = BI   * [D, f];
VtBIDf = BIV' * [D, f];
CtBIDf = C'   * BIDf;
DtBIDf = D'   * BIDf;

M  =     CtBIDf(:, 1:end-1);    % == C' * inv(B) * D
MV = M - VtBIDf(:, 1:end-1);    % == M - V'*inv(B)*D
S  =     DtBIDf(:, 1:end-1);    % == D' * inv(B) * D

g = g - CtBIDf(:, end) ...
      + VtBIDf(:, end);         % == g - C'*inv(B)*f + V'*inv(B)*f
h = h - DtBIDf(:, end);         % == h - D'*inv(B)*f

L       = spdiags(diag(C'*BI*C - BIV'*C - P), 0, ncell, ncell);
LIMVg   = L  \ [MV, g];         % use MLDIVIDE for improved accuracy
MtLIMVg = M' * LIMVg;

S = S - MtLIMVg(:, 1:end-1);    % == S - M.'*inv(L)*MV
h = MtLIMVg(:,end) - h;         % == M.'*inv(L)*g - h

%--------------------------------------------------------------------------
% Solve Schur system to recover unknown contact pressures -----------------
%   (D.'*BI*D - M.'*inv(L)*MV) cp = right hand side
%
lam   = opt.LinSolve(S, h);

%--------------------------------------------------------------------------
% Recover cell pressure from reduced (diagonal) system --------------------
%   L p = g + MV*lam
%
press = LIMVg(:, end) + LIMVg(:, 1:end-1)*lam;

if issparse(press),
   % This happens whenever NUMEL(lam)==1
   press = full(press);
end

%--------------------------------------------------------------------------
% Recover half-contact fluxes from reduced system -------------------------
%   B v = f + C*p - D*lam
%
flux  = BIDf(:, end) - BIDf(:, 1:end-1)*lam + BI*C*press;
