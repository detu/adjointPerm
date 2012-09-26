function [faceFlux, press, lam_n, varargout] = tpf(B, C, D, V, P, f, g, h, ...
                                                         orientation,      ...
                                                         neumannFaces,     ...
                                                         varargin)
%Solve symmetric system of linear eqns using reduction to a mixed system.
%
% SYNOPSIS:
%   [faceFlux, press, lam_n]    = mixedSymm(B, C, D, f, g, h, ...
%                                           orientation, neumannFaces)
%   [faceFlux, press, lam_n]    = mixedSymm(B, C, D, f, g, h, ...
%                                           orientation, neumannFaces, ...
%                                           'pn', pv, ...)
%   [faceFlux, press, lam_n, A] = mixedSymm(...)
%
% DESCRIPTION:
%   Solves the symmetric, hybrid (block) system of linear equations
%
%       [  B     C   D  ] [  v ]     [ f ]
%       [  C'-V' P   0  ] [ -p ]  =  [ g ]
%       [  D'    0   0  ] [ cp ]     [ h ]
%
%   with respect to the face flux `vf', the cell pressure `p', and the face
%   (contact) pressure `cp' for Neumann boundary.  The system is solved
%   using a reduction to a symmetric mixed system of linear equations
%
%       [  Bm   Cm   Dm ] [ vf   ]     [  fm ]
%       [  Cm'  P    0  ] [ -p   ]  =  [  g  ]
%       [  Dm'  0    0  ] [ cp_n ]     [ h_n ]
%
%   where Bm, Cm are 'the' mixed block matrices, and Dm enforces Neumann
%   BCs (size(Dm) = n x nb, where n is the number of faces and nb is the
%   number of boundary faces).
%
% PARAMETERS:
%   B       - Matrix B  in the block system.  Assumed SPD.
%   C       - Block `C' of the block system.
%   D       - Block `D' of the block system.
%   f       - Block `f' of the block system right hand side.
%   g       - Block `g' of the block system right hand side.
%   h       - Block `h' of the block system right hand side.
%   orientation -
%             vector of orientation of half-contacts
%             (i.e. column 2 of G.cellFaces)
%   'pn'/pv - List of name/value pairs describing options to solver.
%             The following options are supported:
%               - 'Regularize' -- TRUE/FALSE, whether or not to enforce
%                                 p(1)==0 (i.e., set pressure zero level).
%                                 This option may be useful when solving
%                                 pure Neumann problems.
%
%               - 'MixedB'     -- TRUE/FALSE (Default FALSE). Must be set
%                                 to TRUE if input 'B' is actually the
%                                 mixed dicsretization matrix Bm (and not
%                                 the hybrid). Useful for multiscale when
%                                 Bm is formed directly.
%
%               - 'Face2CellFace' --
%                                 Sparse matrix which is equal to D by
%                                 default. Differs from D when total rate
%                                 BCs are given.
%
% RETURNS:
%   faceFlux - The face fluxes `vf'.
%   press    - The cell pressure `p'.
%   lam_n    - The contact pressure `cp' on Neumann faces.
%   A        - Mixed system matrix.
%              OPTIONAL.  Only returned if specifically requested.

%{
#COPYRIGHT#
%}

% $Id: mixedSymm.m 529 2008-08-12 09:22:27Z bska $

opt = struct('Regularize',    false, ...
             'MixedB',        false, ...
             'Face2CellFace', []);
opt = merge_options(opt, varargin{:});

regularize = opt.Regularize;
mixedB     = opt.MixedB;
DHat       = opt.Face2CellFace;

numCF = length(orientation);
numC  = size(C, 2);
numNF = nnz(neumannFaces);

% Compute Do-matrix which transforms half-contact fluxes v to contact
% fluxes vf via v = Do * vf

if isempty(DHat)
    Do = spdiags(orientation, 0, numCF, numCF) * D;
    numF = size(D, 2);
else
    Do = spdiags(orientation, 0, numCF, numCF) * DHat;
    numF = size(DHat, 2);
end

if ~mixedB,
   Bm = Do.' * B * Do;
else
   Bm = B;
end
%Cm      = Do.' * C;
%CmHat   = Do.' * (C - V);
%Dm      = Do.' * D(:, neumannFaces);
fm      = Do.' * f;
CmDm    = Do.' * [C, D(:, neumannFaces)];
CmHatDm = Do.' * [C-V, D(:, neumannFaces)];

% Mixed system
% A = [Bm     Cm                    Dm
%      CmHat' P                     sparse(numC,  numNF)
%      Dm'    sparse(numNF, numC)   sparse(numNF, numNF)];
%
% TPF system
invBm = spdiags(1./diag(Bm), 0, numF, numF);
A     = CmHatDm' * invBm * CmDm - spdiags([diag(P); zeros(numNF, 1)], 0, numC + numNF, numC + numNF);
rhs = [g; h(neumannFaces)] - CmHatDm'*invBm*fm;

% Solve system
if regularize
    A(1, 1) = 1;
end
x = A \ rhs;

press     =  x(1:numC);
lam_n     =  x(end-numNF+1 : end);
faceFlux  =  invBm*(CmDm*x + fm);

% faceFlux =  x(              1 : numF       );
% press    = -x(       numF + 1 : numF + numC);
% lam_n    =  x(end - numNF + 1 : end        );

if nargout > 3, varargout{1} = A; end
