function [secAdjRes] = solveSecondOrderAdjointPressureSystem(G, S, W, rock, fluid, simRes, adjRes, linSimRes, secAdjRes, obj, varargin)
cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C     = sparse(1:numel(cellNo), cellNo, 1);

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
curStep = find( cellfun(@(x)~isempty(x), {secAdjRes.timeInterval}), 1, 'first');
dt      = simRes(curStep).timeInterval * [-1 1]';

% Set various
numW     = numel( W );
pvInv    = 1./( G.cells.volumes.*rock.poro );
f_w     = fluid.fw( simRes(curStep).resSol );
dfw      = fluid.dfw(simRes(curStep).resSol);
dltInvm1 = fluid.dLtInv(simRes(curStep-1).resSol);

% Generate RHS, that is f-part (rest is zero)
[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                                        simRes(curStep).wellSol, 'VectorOutput', true);

dQPluss =  double( signQ > 0 );
dQMinus =  double( signQ < 0 );


% Generate right-hand-side ------------------------------------------------
% -------------------------------------------------------------------------

%---------------- First parts similar to standard adjoint -----------------
p_s     = secAdjRes(curStep).resSol.s;
f_bc    =  - dt*( f_w(A.j).*( S.C*(pvInv.*p_s) ) ) ...
           + dt*( S.C*( (f_w.*dQMinus + dQPluss).*( pvInv.*p_s ) ) );


for wellNr = 1 : numel(W)
    W(wellNr).S.RHS.h = zeros( size(W(wellNr).S.RHS.h) );
end

%---------------- Compute d2(F'l)/dx2 w  ----------------------------------
l_v    = adjRes(curStep).resSol.cellFlux;           % l_v^n
l_s    = adjRes(curStep).resSol.s;                 % l_s^n
w_sm1  = linSimRes(curStep-1).resSol.s;            % w_s^{n-1}
w_s    = linSimRes(curStep).resSol.s;              % w_s^n
w_p    = linSimRes(curStep).resSol.s;              % w_p^n -ADDED BY EKA
tpl    = dt * ( pvInv .* l_s); 

% (dsm1, dv) (B(sm1)v, l_v) part (minus - plus)
f_bc = f_bc ...
    - (S.C*dltInvm1) .* (S.B * l_v) .* (S.C*w_sm1);

% (ds, dv) (g(v,s), l_s) part (minus - minus)
f_bc = f_bc ...
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
%--ADDED BY EKA-- Include d2J/dp2 BHP in obj.function----------------------
inx = 0;
for wellNr = 1 : numel(W) %  (plus)
    numCells = length( W(wellNr).cells );
    W(wellNr).S.RHS.f = W(wellNr).S.RHS.f ...
        + obj.partials(curStep).qs(:, inx +(1:numCells))' * w_s;  
%     W(wellNr).S.RHS.h = W(wellNr).S.RHS.h + obj.partials(curStep).p2( W(wellNr).cells ) * w_p( W(wellNr).cells );  % if obj.func. depends on BHP
    inx = inx + numCells;
end


% Solve linear system based on s^{n-1}
b    = computeAdjointRHS(G, W, f_bc);

% % obj. func. is function of pressure in all grid blocks
% b{2} = obj.partials(curStep).p2 * w_p;  % if obj.func. depends on BHP


if strcmp(S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end

%---------------- Solve linear system based on s^{n-1} --------------------

[resSol, wellSol] = solveIncompFlow(simRes(curStep-1).resSol, [], G, ...
                                    S, fluid, 'wells', W, 'rhs', b,  ...
                                    'Solver', solver);

    

%---------------- Solve linear system based on s^{n-1} --------------------
%[resSol, wellSol] = solveWellSystem(simRes(curStep-1).resSol, G, S, W, fluid);

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
