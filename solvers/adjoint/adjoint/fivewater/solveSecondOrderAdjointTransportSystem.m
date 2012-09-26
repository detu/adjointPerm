function [secAdjRes] = solveSecondOrderAdjointTransportSystem(G, S, W, rock, fluid, simRes, adjRes, linSimRes, secAdjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
numSteps = numel(simRes);
if isempty(secAdjRes)
    curStep = numSteps;
else
    curStep = find( cellfun('isempty', {secAdjRes.timeInterval}), 1, 'last');
end
dt      = simRes(curStep).timeInterval * [-1 1]';
secAdjRes(curStep).timeInterval = simRes(curStep).timeInterval;


% Set various
numC     = G.cells.num;
numCF    = size(S.B, 1);
numW     = numel( W );
pvInv    = 1./( G.cells.volumes.*rock.poros );
PvInv    = spdiags(pvInv, 0, numC, numC);

dfw      = fluid.Dfw(simRes(curStep).resSol.sw);
Dfw      = spdiags(dfw, 0, numC, numC);
d2fw     = fluid.D2fw(simRes(curStep).resSol.sw);
dltInv   = fluid.DLtInv(simRes(curStep).resSol.sw);
d2ltInv  = fluid.D2LtInv(simRes(curStep).resSol.sw);

%--------------------------------------------------------------------------

% Generate system matrix
[A, qPluss, signQ]  = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                              simRes(curStep).wellSol, 'VectorOutput', true);

At      = sparse(A.j, A.i, -simRes(curStep).resSol.cellFlux, numC, numC) ...
            + spdiags(A.qMinus, 0, numC, numC);

systMat = speye(numC, numC) - dt * ( Dfw * At * PvInv);   % system matrix

% Some extra A - related
dQPluss =  double( signQ > 0 );
dQMinus =  double( signQ < 0 );
Cj      = sparse( (1:numCF)', A.j, 1, numCF, numC);

% Generate right-hand-side ------------------------------------------------
% -------------------------------------------------------------------------
RHS = zeros(numC, 1);

%---------------- First parts similar to standard adjoint -------------
if (curStep < numSteps)
    v1       = simRes(curStep+1).resSol.cellFlux;       % v^{n+1}
    p_v1     = secAdjRes(curStep+1).resSol.cellFlux;    % p_v^{n+1}
    
    RHS     = secAdjRes(curStep+1).resSol.sw;
    
    % d/ds B^{n+1}v^{n+1} part
    RHS     = RHS - S.C'*spdiags( (S.C*dltInv).*v1 , 0, numCF, numCF)*S.B*p_v1;
    
    % d/ds B_w^{n+1}q_w^{n+1} part
    for wellNr = 1:numW
        w   = W(wellNr);
        q1  = simRes(curStep+1).wellSol(wellNr).flux;       % q_w^{n+1}
        p_q1= secAdjRes(curStep+1).wellSol(wellNr).flux;    % p_q_w^{n+1}
        RHS = RHS + w.S.C'*( ( (w.S.C*dltInv).*q1 ).*(w.S.B*p_q1) );
    end
end

%---------------- Compute ddFl =d2(F'l)/dx2 w  ----------------------------
ddFl = zeros(numC, 1);
l_s = adjRes(curStep).resSol.sw;                    % l_s^{n}
w_s = linSimRes(curStep).resSol.sw;               % w_s^{n}
w_v = linSimRes(curStep).resSol.cellFlux;           % w_v^{n}
tpl = dt * ( pvInv .* l_s); 

if (curStep < numSteps)  
    l_v1 = adjRes(curStep+1).resSol.cellFlux;            % l_v^{n+1}
    w_v1 = linSimRes(curStep+1).resSol.cellFlux;         % w_v^{n+1}
    
    % (dv1, ds) (B(s)v1, l_v1) part (plus)
    ddFl = ddFl + S.C'*( (S.C*dltInv) .* l_v1 .*(S.B*w_v1) );
     
    % (ds, ds) (B(s)v1, l_v1) part (plus)
    ddFl = ddFl + S.C' * ( (S.C*d2ltInv) .* l_v1 .* (S.B*v1) .* (S.C*w_s) );
 
    for wellNr = 1:numW
        w    = W(wellNr);
        q1   = simRes(curStep+1).wellSol(wellNr).flux;     % q_w^{n+1}
        w_q1 = linSimRes(curStep+1).wellSol(wellNr).flux;  % w_{q_w}^{n+1}
        l_q1 = adjRes(curStep+1).wellSol(wellNr).flux;     % l_{q_w}^{n+1}
        
        % (dq1, ds) (B_w(s)q1, l_q1) part (minus)
        ddFl = ddFl - w.S.C'*( (w.S.C*dltInv) .* l_q1 .*(w.S.B*w_q1) );
        
        % (ds, ds) (B_w(s)q1, l_q1) part (minus)
        ddFl = ddFl - w.S.C' ...
                * ( (w.S.C*d2ltInv) .* l_q1 .* (w.S.B*q1) .* (w.S.C*w_s) );
     end
end
   
% (dv, ds) (g(v,s), l_s) part (minus)
ddFl = ddFl - (- Cj' * ( dfw(A.j) .* (S.C*tpl) .* w_v) ...
        + spdiags( dfw.* dQMinus .* tpl, 0, numC, numC) * (S.C'*w_v) );
 
% (ds, ds) (g(v,s), l_s) part (minus)
ddFl = ddFl - d2fw .* (At*tpl) .* w_s;

% Finally update RHS:
RHS = RHS - ddFl;

%----------  d2/ds2 and d2/dsdq of objective ------------------------------
ddo = obj.partials(curStep).s2 * w_s; 

inx = 0;
for wellNr = 1 : numel(W)
    w_q1 = linSimRes(curStep).wellSol(wellNr).flux;  % w_{q_w}^n
    numCells = length( W(wellNr).cells );
    ddo = ddo + obj.partials(curStep).qs(:, inx + (1:numCells)) * w_q1;
    inx = inx + numCells;
end

RHS = RHS - ddo; 

%---------- Solve system --------------------------------------------------
secAdjRes(curStep).resSol.sw  =  systMat \ RHS;

return
