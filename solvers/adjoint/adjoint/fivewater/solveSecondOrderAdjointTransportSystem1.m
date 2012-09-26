function [secAdjRes] = solveSecondOrderAdjointTransportSystem1(G, S, W, rock, fluid, simRes, adjRes, linSimRes, secAdjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
numSteps = numel(simRes);
if isempty(secAdjRes)
    curStep = numSteps;
else
    curStep = find( cellfun('isempty', {secAdjRes.timeInterval}), 1, 'last');
end
dt      = simRes(curStep).timeInterval * [-1 1]';
secAdjRes(curStep).timeInterval = simRes(curStep).timeInterval;

% Generate system matrix
numC    = G.cells.num;
numCF   = size(S.B, 1);
PV      = G.cells.volumes.*rock.poros;
invDPV  = spdiags(1./PV, 0, numC, numC);
DDf     = spdiags( fluid.Dfw(simRes(curStep).resSol.sw), 0, numC, numC);
DD2f    = spdiags( fluid.D2fw(simRes(curStep).resSol.sw), 0, numC, numC);
DDLtI   = spdiags( fluid.DLtInv(simRes(curStep).resSol.sw), 0, numC, numC);
DD2LtI  = spdiags( fluid.D2LtInv(simRes(curStep).resSol.sw), 0, numC, numC);

[A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                              simRes(curStep).wellSol, 'VectorOutput', true);

At      = sparse(A.j, A.i, -simRes(curStep).resSol.cellFlux, numC, numC) + spdiags(A.qMinus, 0, numC, numC);

systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
%clear PV invDPV DDf At

% Generate right-hand-side

% d2/ds2 of objective
w_s_n = linSimRes(curStep).resSol.sw;      % w_s^{n}
RHS   = -obj.partials(curStep).s2 * w_s_n; 

numW    = numel(W);
dim     = fluid.DLtInv(simRes(curStep).resSol.sw);
d2im    = fluid.D2LtInv(simRes(curStep).resSol.sw);

l_s     = adjRes(curStep).resSol.sw;   % l_s^{n}
    
if (curStep < numSteps)
    w_v_n1  = linSimRes(curStep+1).resSol.cellFlux;      % w_v^{n+1}
    l_v_n1  = adjRes(curStep+1).resSol.cellFlux;   % l_v^{n+1}
    v       = simRes(curStep+1).resSol.cellFlux;   % v^{n+1}
    p_v     = secAdjRes(curStep+1).resSol.cellFlux;   % p_v^{n+1}
   
    RHS     = secAdjRes(curStep+1).resSol.sw;
    
    % Update d/ds B^{n+1}v^{n+1} part
    RHS     = RHS - S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.B*p_v;
    
    % Update d2/ds2 B^{n+1}v^{n+1} part
    RHS = RHS - d2im .* ( S.C'*spdiags(v,0,numCF,numCF)*S.B*l_v_n1 ) .* w_s_n;
    
    % Update d2/dv{n+1}ds B^{n+1}v^{n+1} part
    RHS = RHS - S.C'* spdiags(l_v_n1, 0, numCF, numCF) * S.B * ( (S.C*dim).*w_v_n1 );   
   
    % Update d/ds B_w^{n+1}q_w^{n+1} part
    for wellNr = 1:numW
        w   = W(wellNr);
        q   = simRes(curStep+1).wellSol(wellNr).flux; % q_w^{n+1}
        p_q = secAdjRes(curStep+1).wellSol(wellNr).flux; % p_q_w^{n+1}
        RHS = RHS + w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.B*p_q) );
    end
    
    % Update q_w parts
    for wellNr = 1:numW
        w   = W(wellNr);
        w_qw_n1 = linSimRes(curStep+1).wellSol(wellNr).flux;
        l_qw_n1 = adjRes(curStep+1).wellSol(wellNr).flux;
        RHS = RHS - w.S.C'* w.S.B * ( (w.S.C*dim).*l_qw_n1.*w_qw_n1 );
        
        q   = simRes(curStep+1).wellSol(wellNr).flux; % q_w^{n+1}
%         RHS = RHS + d2im .* (w.S.C' * sparse(diag(q)) * w.S.B * l_qw_n1) .* w_s_n;
        RHS = RHS + d2im .* (w.S.C' * q * w.S.B * l_qw_n1) .* w_s_n;
    end
end
   

%---------- d2/ds2 for res eq
% Update d2/ds2 B^{n+1}v^{n+1} part

% RHS = RHS - d2im .* ( S.C'*spdiags(v,0,numCF,numCF)*S.B*l_v_n1 ) .* w_s_n;

% Update d2g/ds2 (g = nonlinear part of saturation eq)
% RHS = RHS + dt * DD2f * ( (l_s'* invDPV * At)' * w_s_n);

RHS = RHS + dt * invDPV * DD2f * ( ( l_s' *At)' .* w_s_n);

% % Update q_w parts
% for wellNr = 1:numW
%     w   = W(wellNr);
% %     l_qw_n1 = adjRes(curStep+1).wellSol(wellNr).flux;
% %     q   = simRes(curStep+1).wellSol(wellNr).flux; % q_w^{n+1}
%     
%     if curStep == numSteps
%         l_qw_n1 = zeros( size( adjRes(curStep).wellSol(wellNr).flux) );
%         q       = zeros( size( simRes(curStep).wellSol(wellNr).flux) ); % q_w^{n+1}
%     else
%         l_qw_n1 = adjRes(curStep+1).wellSol(wellNr).flux;
%         q       = simRes(curStep+1).wellSol(wellNr).flux; % q_w^{n+1}
%     end
%             
%     
%     RHS = RHS + d2im .* (w.S.C' * sparse(diag(q)) * w.S.B * l_qw_n1) .* w_s_n;
% end

% -----------------

% Update d2/dvds
Cu = sparse( (1:numCF)', A.j, 1, numCF, numC );
dQPluss =  double( signQ > 0 );
dQMinus =  double( signQ < 0 );
DDQP    =  spdiags(dQPluss, 0, numC, numC);
DDQM    =  spdiags(dQMinus, 0, numC, numC);
l_s     =  adjRes(curStep).resSol.sw;
w_v_n   =  linSimRes(curStep).resSol.cellFlux;      % w_v^{n}

RHS  = RHS - dt * spdiags(l_s, 0, numC ,numC) * invDPV * DDf * Cu' * w_v_n ...
           + dt * spdiags(l_s, 0, numC ,numC) * invDPV * ( DDf*DDQP + DDQM) * S.C' * w_v_n;



% Solve system
secAdjRes(curStep).resSol.sw  = systMat \ RHS;

return
