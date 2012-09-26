function [adjRes] = solveAdjointPressureSystemMS(G, CG, S, CS, W, rock, fluid, simRes, adjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
curStep = find( cellfun(@(x)~isempty(x), {adjRes.timeInterval}), 1, 'first');
dt      = simRes(curStep).timeInterval * [-1 1]';

% get partition
 p     = zeros([size(CG.cells.subCells,1), 1]);
 [i,j] = find(CG.cells.subCells);
 p(i)  = j;
      
wellSol = initWellSol(W, 0);

% Generate RHS, that is f-part (rest is zero)
PV      = G.cells.volumes.*rock.poro;
invPV   = 1./PV;
%v       = simRes(curStep).resSol.cellFlux;  % v^n
%q       = S.C'*v;                           % q^n
%s_qp    = (q > 1e-10);                      % sign(q^n_+)  XXXXXXXX Make threshold relative XXXXXX              
%s_qm    = -(q < -1e-10);                    % sign(q^n_-)  XXXXXXXx
%f_w     = fluid.krw(simRes(curStep).resSol) ./ fluid.Lt(simRes(curStep).resSol); % fractinal flow, f(s^n)
f_w     = fluid.fw( simRes(curStep).resSol);
l_s     = adjRes(curStep).resSol.s;        % lam_s^n
% Flux-matrix: A.i, A.j, A.qMinus
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                                        simRes(curStep).wellSol, 'VectorOutput', true);
cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);

dQPluss =  double( signQ > 0 );
dQMinus = -double( signQ < 0 );

f_bc    = -obj.partials(curStep).v';
f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
               + dt*( S.C*( (-f_w.*dQMinus + dQPluss).*( invPV.*l_s ) ) );
%f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
%               + dt*( S.C*( (f_w.*s_qm + s_qp).*( invPV.*l_s ) ) );


% Assume all f, g, h are zero  in S.RHS, set f_bc as above. Set f_bc = 0
% after solve
% S.RHS.f_bc  = f_bc;
% if strcmp(S.type, 'mixed')
%     cellNo      = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
%     orientation = 2*(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo)-1;
%     %no longer works: orientation = G.cellFaces(:,2); 
%     numcf       = numel(orientation);
%     Do          = spdiags(orientation, 0, numcf, numcf) * S.D;
%     r_fac       = Do'*f_bc;
% else
%     r_fac = S.BI*f_bc;
% end

mob   = fluid.Lt(simRes(curStep).resSol);
omega = fluid.omega(simRes(curStep).resSol);

[B, C, D, f, h, Psi, Phi] = ...
      unpackWellSystemComponentsMS(W, G, p, mob, omega);

Psi = S.BI*Psi;     

% Set f and h in W(i).S.RHS equal DJ/Dp_w and zero. Leave naumannFaces and dirichletFaces as is 
inx = 0;
for wellNr = 1 : numel(W)
    numCells = length( W(wellNr).cells );
    W(wellNr).S.RHS.f  = obj.partials(curStep).q_w( inx+1 : inx+numCells )';
    %W(wellNr).CS.RHS.f = W(wellNr).CS.rates'*W(wellNr).S.RHS.f - W(wellNr).CS.basis'*r_fac;
    W(wellNr).CS.RHS.f = W(wellNr).CS.rates'*W(wellNr).S.RHS.f - Psi(:,wellNr)'*f_bc;
    W(wellNr).S.RHS.h  = zeros( size(W(wellNr).S.RHS.h) );
    W(wellNr).CS.RHS.h = zeros( size(W(wellNr).CS.RHS.h) );
    inx = inx + numCells;
end

if strcmp(S.type, 'hybrid')
   solver = 'hybrid';   
else
   solver = 'mixed';
end

b = computeAdjointRHSMS(G, CG, S, CS, W, f_bc);


% Solve linear system based on s^{n-1}
[resSol, wellSol] = solveIncompFlowMS(simRes(curStep-1).resSol, wellSol, ...
                                      G, CG, p, S, CS, fluid, 'wells', W,...
                                      'rhs', b, 'Solver', solver);


% Update adjRes !!! Note minuses in front of pressure and wellrates in
% forward system, but not in adjoint, thus set minus here  
adjRes(curStep).resSol.cellFlux     = resSol.cellFlux;
adjRes(curStep).resSol.cellPressure = - resSol.cellPressure;           % !!!minus
adjRes(curStep).resSol.facePressure = resSol.facePressure;
adjRes(curStep).wellSol             = wellSol;
for k = 1 : numel(adjRes(curStep).wellSol)
    adjRes(curStep).wellSol(k).flux = - adjRes(curStep).wellSol(k).flux;     % !!!minus
end
end
%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------
 
function b = computeAdjointRHSMS(G, CG, S, CS, W, f_res) 
% Compute adjoint 'pressure' rhs 
% for use in function solveAdjointPressureSystem
%
% SYNOPSIS:
%   b = computeAdjoint(G, W, f_res, f_w)
% 
% DESCRIPTION:
% Computes adjoint 'pressure' rhs as input to solveIncompFlow 
% in function solveAdjointPressureSystem
%
% PARAMETERS:
%   G     - Grid data structure.
%  
%   W     - Well structure as defined by addWell &c.
%
%   f_res - Adjoint reservoir 'pressure' condtions
%
% RETURNS:
%   b     - Ajoint pressure rhs to be passed directly as option
%           'rhs' to solveIncompFlow.
%

assert( numel(f_res) == size(G.cellFaces,1)); 
f_w = [];

if strcmpi(S.type, 'hybrid'),
   [Bv, Phi] = basisMatrixHybrid(G, CG, CS);
else
   [Bv, Phi] = basisMatrixMixed (G, CG, CS); 
end

Psi = S.BI * Bv;
f_res = Psi.' * f_res;

% unpack adjoint well 'pressure' conditions
if ~isempty(W),  
   CS_w  = [ W.CS   ];
   RHS = [ CS_w.RHS ];
   f_w = {RHS.f};
   f_w = vertcat(f_w{:});
   h_w = {RHS.h};
   h_w = vertcat(h_w{:});    
end

% b = [f_res; f_w; g; h_res; h_w]
b    = cell([1, 3]);
b{1} = vertcat(f_res, f_w);
b{2} =         zeros([CG.cells.num, 1]);
b{3} = vertcat(zeros([size(CS.activeFaces,1), 1]), h_w);
end


