function [adjRes] = solveAdjointTransportSystem(G, S, W, rock, fluid, simRes, adjRes, obj, varargin)

% Find current time step (search for empty slots in adjRes)
% NOTE: actually curent time step +1
numSteps = numel(simRes);
if isempty(adjRes)
    curStep = numSteps;
else
    curStep = find( cellfun('isempty', {adjRes.timeInterval}), 1, 'last');
end
dt      = simRes(curStep).timeInterval * [-1 1]';
adjRes(curStep).timeInterval = simRes(curStep).timeInterval;

% Generate system matrix
numC    = G.cells.num;
PV      = G.cells.volumes.*rock.poro;
invDPV  = spdiags(1./PV, 0, numC, numC);
DDf     = spdiags( fluid.dfw(simRes(curStep).resSol), 0, numC, numC);
At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                              simRes(curStep).wellSol, 'Transpose', true);
systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
clear PV invDPV DDf At

% Generate right-hand-side
RHS     = -obj.partials(curStep).s';
if curStep < numSteps
    numCF   = size(S.B, 1);
    numW    = numel(W);
    RHS     = RHS + adjRes(curStep+1).resSol.s;
    dim     = fluid.dLtInv(simRes(curStep).resSol);
    
    % Update B^{n+1}v^{n+1} part
    v       = simRes(curStep+1).resSol.cellFlux;   % v^{n+1}
    l_v     = adjRes(curStep+1).resSol.cellFlux;   % lam_v^{n+1}
    cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
    S.C     = sparse(1:numel(cellNo), cellNo, 1);
    RHS     = RHS - S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.B*l_v;
    
    % Update B_w^{n+1}q_w^{n+1} part
    for wellNr = 1:numW
        w   = W(wellNr);
        q   = simRes(curStep+1).wellSol(wellNr).flux; % q_w^{n+1}
        l_q = adjRes(curStep+1).wellSol(wellNr).flux; % lam_q_w^{n+1}
        RHS = RHS + w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.B*l_q) );
    end
end    

% Solve system
adjRes(curStep).resSol.s  = systMat \ RHS;

return

