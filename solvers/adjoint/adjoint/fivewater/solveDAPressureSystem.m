function [adjRes] = solveDAPressureSystem(G, S, W, rock, fluid, simRes, adjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
curStep = find( cellfun(@(x)~isempty(x), {adjRes.timeInterval}), 1, 'first');
dt      = simRes(curStep).timeInterval * [-1 1]';

% Generate RHS, that is f-part (rest is zero)
PV      = G.cells.volumes.*rock.poros;
invPV   = 1./PV;
%f_w     = fluid.krw(simRes(curStep).resSol) ./ fluid.Lt(simRes(curStep).resSol); % fractinal flow, f(s^n)
f_w     = fluid.fw( simRes(curStep).resSol.sw );
l_s     = 0*adjRes(curStep).resSol.sw;        % lam_s^n
% Flux-matrix: A.i, A.j, A.qMinus
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                                        simRes(curStep).wellSol, 'VectorOutput', true);

dQPluss =  double( signQ > 0 );
dQMinus = -double( signQ < 0 );

f_bc    = -obj.partials(curStep).v';
f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
               + dt*( S.C*( (-f_w.*dQMinus + dQPluss).*( invPV.*l_s ) ) );
%f_bc    = f_bc - dt*( f_w(A.j).*( S.C*(invPV.*l_s) ) ) ...
%               + dt*( S.C*( (f_w.*s_qm + s_qp).*( invPV.*l_s ) ) );

S.RHS.f_bc  = f_bc;

% Set f and h in W(i).S.RHS equal DJ/Dp_w and zero. Leave naumannFaces and dirichletFaces as is 
inx = 0;
for wellNr = 1 : numel(W)
    numCells = length( W(wellNr).cells );
    W(wellNr).S.RHS.f = obj.partials(curStep).q_w( inx+1 : inx+numCells )';
    W(wellNr).S.RHS.h = zeros( size(W(wellNr).S.RHS.h) );
    inx = inx + numCells;
end

% Solve linear system based on s^{n-1}
[resSol, wellSol] = solveWellSystem(simRes(curStep-1).resSol, G, S, W, fluid);

% Update adjRes !!! Note minuses in front of pressure and wellrates in
% forward system, but not in adjoint, thus set minus here  
adjRes(curStep).resSol.cellFlux     = resSol.cellFlux;
adjRes(curStep).resSol.cellPressure = - resSol.cellPressure;           % !!!minus
adjRes(curStep).resSol.facePressure = resSol.facePressure;
adjRes(curStep).wellSol             = wellSol;
for k = 1 : numel(adjRes(curStep).wellSol)
    adjRes(curStep).wellSol(k).flux = - adjRes(curStep).wellSol(k).flux;     % !!!minus
end
    

return;
