function [resSol, wellSol, S] = solveTPFBlackOilWellSystem(resSol, wellSol, ...
                                                        G, rock, S, W,   ...
                                                        PVTTAB, p0, dt,  ...
                                                        varargin)
% solveTPFBlackOilWellSystem -- Perform one iteration of successive substitution
%                               algorithm on Black Oil pressure system (including
%                               wells). This is just hacked directly from 
%                               solveBlackOilWellSystem, and does not use
%                               mixed type system as one would expect.
%
% SYNOPSIS:
%   [resSol, wellSol] = solveBlackOilWellSystem(resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt)
%   [resSol, wellSol] = solveBlackOilWellSystem(resSol, wellSol, ...
%                                               G, rock, S, W,   ...
%                                               PVTTAB, p0, dt,  ...
%                                               'pn1', pv1, ...)
%
% PARAMETERS:
%   resSol  - Reservoir solution structure from previous time step (or
%             previous iteration of successive substitution algorithm).
%
%   wellSol - Well solution structure from previous time step (or previous
%             iteration of successive substitution algorithm).
%
%   G       - Grid structure as described in grid_structure.
%
%   rock    - Rock structure.
%
%   S       - Linear system structure as defined by function
%             assembleMixedSystem.
%
%   W       - Well linear system structure as defined by function
%             assembleWellSystem.
%
%   PVTTAB  - PVT tables.
%
%   p0      - Vector, length G.cells.num, of cell pressures at previous time
%             step.
%
%   dt      - Time step size.
%
%   'pn'/pv - List of name/value pairs encompassing the union of optional
%             parameters to functions updateBlackOilMimeticSystem and
%             updateBlackOilWellSystem.
%
% RETURNS:
%   resSol  - Updated reservoir solution structure.
%
%   wellSol - Updated well solution structure.

% $Id: solveBlackOilWellSystem.m 771 2008-09-23 17:16:27Z bska $

%--------------------------------------------------------------------------
% Get updated pressure and sat-dependent quantities -----------------------
%

[fluid, resSol] = pvt(PVTTAB, resSol, S.constants.gravity);

S = updateBlackOilMimeticSystem(G, S, rock, fluid, p0, resSol, dt, ...
                                varargin{:});
W = updateBlackOilWellSystem   (G, S, W,    fluid, wellSol, ...
                                varargin{:});

%--------------------------------------------------------------------------
% Extract well system contributions ---------------------------------------
%
[BIW, BIVW, CW, DW, fW, hW] = unpackWellSystemComponents(W, fluid.Lt);

%--------------------------------------------------------------------------
% Build hybrid system components ------------------------------------------
%

DMob = spdiags(S.C * fluid.Lt, 0, S.sizeB(1), S.sizeB(2));
BI   = blkdiag(DMob * S.BI , BIW {:});
BIV  = vertcat(       S.BIV, BIVW{:});
C    = vertcat(       S.C  , CW  {:});
D    = blkdiag(       S.D  , DW  {:});
P    = S.P;

DHat = blkdiag(S.D, diag(vertcat(DW{:})));

% ---- Hack -----------------------
B    = spdiags(1./diag(BI), 0, size(BI, 1), size(BI, 1));
V    = B*BIV;
clear BI BIV
% ---------------------------------


f_mimetic = S.RHS.f_bc + ...
            S.RHS.f_grav .* (S.C * fluid.omega) .* S.constants.gravity;

f    = vertcat(f_mimetic, fW{:});
g    = S.RHS.g_src + S.RHS.g_comp + S.RHS.g_grav + ...
       S.RHS.volume_discrepancy;
h    = vertcat(S.RHS.h_bc, hW{:});

%--------------------------------------------------------------------------
% Solve global hybrid system (unsymmetric) --------------------------------
%
orientation  = [G.cellFaces(:,2); ...
                ones([size(B,1) - S.sizeB(1), 1])];
           
neuFW  = strcmp({ W.type } .', 'rate');
neuFR  = S.RHS.neumannFaces;
numNFR = nnz(neuFR);
neuF   = [neuFR; neuFW];

% dirFR    =  S.RHS.dirichletFaces;
% dirCondR = -S.RHS.f_bc(logical(S.D * dirFR));
% if ~isempty(W),
%    dirFW    = strcmp({ W.type } .', 'bhp');
%    dirCondW = [ W(dirFW).val ].';
% else
%    dirFW    = logical([]);
%    dirCondW = [];
% end
% dirF    = [dirFR;    dirFW];
% dirCond = [dirCondR; dirCondW];
% 
% lam       = zeros(size(dirF));
% lam(dirF) = dirCond;
% [flux, p, lam(~dirF)] = schurComplement(BI, C, D(:,~dirF), BIV, P, ...
%                                          f, g, h(  ~dirF));

[flux, p, lam_n] = tpf(B, C, D, V, P, f, g, h,  ...
                       orientation, neuF, ...
                      'Face2CellFace', DHat);
%--------------------------------------------------------------------------
% Package solution in accessible form -------------------------------------
%
resSol.cellPressure(:)     = p;
%resSol.facePressure(neuFR) = - lam_n(1 : numNFR);
resSol.cellFlux(:)         = faceFlux2cellFlux(G, flux(1 : S.sizeD(2)));
resSol.faceFlux(:)         = flux(1 : S.sizeD(2));

% Obtain face pressures :
pseudoInvD      = spdiags(1./sum(S.D)', 0, size(S.D, 2), size(S.D, 2))*S.D';
numRF = S.sizeD(1);
resSol.facePressure   = pseudoInvD * (f_mimetic - B(1:numRF, 1:numRF)*resSol.cellFlux + S.C*p);

lamW( neuFW) = -lam_n(numNFR + 1 : end);
lamW(~neuFW) = [ W(~neuFW).val ] .';

wellSol = packageWellSol(flux(S.sizeD(2) + 1 : end), lamW, fW, hW);


% resSol.cellPressure(:) = p;
% resSol.facePressure(:) = lam (1 : S.sizeD(2));
% resSol.cellFlux(:)     = flux(1 : S.sizeB(1));
% resSol.faceFlux(:)     = cellFlux2faceFlux(G, flux(1 : S.sizeB(1)));
% 
% lamW    = lam(S.sizeD(2) + 1 : end);
% wellSol = packageWellSol(flux(S.sizeB(1) + 1 : end), lamW, fW, hW);
