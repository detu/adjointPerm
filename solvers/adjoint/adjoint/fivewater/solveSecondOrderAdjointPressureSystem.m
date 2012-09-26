function [secAdjRes] = solveSecondOrderAdjointPressureSystem(G, S, W, rock, fluid, simRes, adjRes, linSimRes, secAdjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
curStep = find( cellfun(@(x)~isempty(x), {secAdjRes.timeInterval}), 1, 'first');
dt      = simRes(curStep).timeInterval * [-1 1]';

% Set various
numW     = numel( W );
pvInv    = 1./( G.cells.volumes.*rock.poros );
f_w     = fluid.fw( simRes(curStep).resSol.sw );
dfw      = fluid.Dfw(simRes(curStep).resSol.sw);
dltInvm1 = fluid.DLtInv(simRes(curStep-1).resSol.sw);

% Generate RHS, that is f-part (rest is zero)
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                                        simRes(curStep).wellSol, 'VectorOutput', true);

dQPluss =  double( signQ > 0 );
dQMinus =  double( signQ < 0 );


% Generate right-hand-side ------------------------------------------------
% -------------------------------------------------------------------------

%---------------- First parts similar to standard adjoint -----------------
p_s     = secAdjRes(curStep).resSol.sw;
f_bc    =  - dt*( f_w(A.j).*( S.C*(pvInv.*p_s) ) ) ...
           + dt*( S.C*( (f_w.*dQMinus + dQPluss).*( pvInv.*p_s ) ) );

S.RHS.f_bc  = f_bc;

for wellNr = 1 : numel(W)
    W(wellNr).S.RHS.h = zeros( size(W(wellNr).S.RHS.h) );
end

%---------------- Compute d2(F'l)/dx2 w  ----------------------------------
l_v    = adjRes(curStep).resSol.cellFlux;           % l_v^n
l_s    = adjRes(curStep).resSol.sw;                 % l_s^n
w_sm1  = linSimRes(curStep-1).resSol.sw;            % w_s^{n-1}
w_s    = linSimRes(curStep).resSol.sw;              % w_s^n
tpl    = dt * ( pvInv .* l_s); 

% (dsm1, dv) (B(sm1)v, l_v) part (minus - plus)
S.RHS.f_bc = S.RHS.f_bc ...
    - (S.C*dltInvm1) .* (S.B * l_v) .* (S.C*w_sm1);

% (ds, dv) (g(v,s), l_s) part (minus - minus)
S.RHS.f_bc = S.RHS.f_bc ...
      - dfw(A.j) .* (S.C*tpl) .* w_s(A.j) ...
        + S.C * ( dfw .* dQMinus .* tpl .* w_s) ;

% (dsm1, dq) (B_w(sm1)q_w, l_q) part (plus - minus) 
for wellNr = 1:numW
    w    = W(wellNr);
    l_q  = adjRes(curStep).wellSol(wellNr).flux;     % l_{q_w}^n
    W(wellNr).S.RHS.f = ...
        - (w.S.C*dltInvm1) .* (w.S.B * l_q) .* (w.S.C*w_sm1) ;
end

%---------------- Include d2/dq_wds of objective---------------------------
inx = 0;
for wellNr = 1 : numel(W) %  (plus)
    numCells = length( W(wellNr).cells );
    W(wellNr).S.RHS.f = W(wellNr).S.RHS.f ...
        + obj.partials(curStep).qs(:, inx +(1:numCells))' * w_s;  
    inx = inx + numCells;
end

%---------------- Solve linear system based on s^{n-1} --------------------
[resSol, wellSol] = solveWellSystem(simRes(curStep-1).resSol, G, S, W, fluid);

% Update adjRes !!! Note minuses in front of pressure and wellrates in
% forward system, but not in adjoint, thus set minus here  
secAdjRes(curStep).resSol.cellFlux     = resSol.cellFlux;
secAdjRes(curStep).resSol.cellPressure = - resSol.cellPressure;           % !!!minus
secAdjRes(curStep).resSol.facePressure = resSol.facePressure;
secAdjRes(curStep).wellSol             = wellSol;
for k = 1 : numel(secAdjRes(curStep).wellSol)
    secAdjRes(curStep).wellSol(k).flux =   - secAdjRes(curStep).wellSol(k).flux;     % !!!minus
end
    

return;
