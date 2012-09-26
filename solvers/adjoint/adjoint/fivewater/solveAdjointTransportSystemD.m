function [adjRes] = solveAdjointTransportSystemD(G, S, W, rock, fluid, simRes, adjRes, obj, simResDer, varargin)

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
PV      = G.cells.volumes.*rock.poros;
invDPV  = spdiags(1./PV, 0, numC, numC);
DDf     = spdiags( fluid.Dfw(simRes(curStep).resSol.sw), 0, numC, numC);
At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep).resSol, ...
                              simRes(curStep).wellSol, 'Transpose', true);
systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
clear PV invDPV DDf At

% Generate right-hand-side
RHS     = -simResDer(curStep).resSol.sw'*obj.partials(curStep).s2';
if curStep < numSteps
    numCF   = size(S.B, 1);
    numW    = numel(W);
    RHS     = RHS + adjRes(curStep+1).resSol.sw;
    dim     = fluid.DLtInv(simRes(curStep).resSol.sw);
    
    % Update B^{n+1}v^{n+1} part
    v       = simRes(curStep+1).resSol.cellFlux;   % v^{n+1}
    l_v     = adjRes(curStep+1).resSol.cellFlux;   % lam_v^{n+1}
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
adjRes(curStep).resSol.sw  = systMat \ RHS;

return

% Functions below not used anymore!!

function [df] = DFracFlow(s, fluid)
% Derivative of fractional flow function of the form
%
%            s²/µw              s²                µw
%    f(s) = ---------------- = ------------ ,  mr=---
%            s²/µw+(1-s)²/µo    s²+mr*(1-s)²      µo
%
%
%            2 mr*s(1-s)
%   df(s) = ---------------
%           (s²+mr*(1-s)²)²        

mr  = fluid.muw/fluid.muo; % !!!!!!!!!!!!!!!!!
df  = ( (2*mr) * s .* (1-s) ) ./ ( (s.^2 + mr*(1-s).^2).^2 );
return

function [dim] = DInvMob(s, fluid)
% Derivative of inverse mobility function
%
%            
%    mob(s)     = s²/µw+(1-s)²/µo    
%
%         

muw  = fluid.muw;
muo  = fluid.muo;
dim  = ( (-2/muw)*s + (2/muo)*(1-s) )./( ( (1/muw)*s.^2 + (1/muo)*(1-s).^2 ).^2 );
return