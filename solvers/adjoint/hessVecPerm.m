function [Hv] = hessVecPerm(G, S, W, rock, fluid, simRes, adjRes, schedule, controls, objectiveFunction, param, varargin)
% hessVec -- Compute Hessian times vector
%
% SYNOPSIS:
%   Hv = hessVec(G, S, W, rock, fluid, simRes, adjRes, schedule, controls, objectiveFunction, pn, pv ...)
%
% DESCRIPTION:
%   Compute hessian times a vector by forward linearized eqs and 2.order
%   adjoint. If vector is not given in varargin (or 'Varargin' field is set to empty),
%   the vector is assumed to be incorporatred in schedule
%   
% PARAMETERS:
%
%
% RETURNS:
%   Hv      - vector of length equal to #controls
%
%
% SEE ALSO:
%  
opt     = struct('Verbose',  false , ...
                 'Vector', []);
opt     = merge_options(opt, varargin{:});
if opt.Verbose
    verboseLevel = 1; 
else
    verboseLevel = 0;
end

% if ~isempty(opt.Vector)
%     dispif(verboseLevel==1, 'Uptating schedule to incorporate vector'); 
%     controls  = updateControls(controls, opt.Vector);
%     schedule  = setToZero(schedule);
%     schedule  = updateSchedule(controls, schedule);
% end

% linSimRes = runLinearized(simRes, G, S, W, rock, fluid, schedule, 'VerboseLevel', verboseLevel);
deltaU    = opt.Vector;
linSimRes = runLinearizedPerm1(simRes, G, S, W, rock, fluid, schedule, deltaU,  'VerboseLevel', 0);

secAdjRes = runSecondOrderAdjointPerm(simRes, adjRes, linSimRes, G, S, W, rock, fluid, ...
                                  schedule, controls, objectiveFunction, param, 'VerboseLevel', verboseLevel);

Hv  = compHvPerm(param, G, S, W, rock, simRes, secAdjRes, schedule, controls, fluid, objectiveFunction);                              
return
                              
%--------------------------------------------------------------------------

function Hv = compHvPerm(param, G, S, W, rock, simRes, secAdjRes, schedule, controls, fluid, objectiveFunction)
% contribution from the obj.func:  partials(step).u2
% from the model: l_v and l_qw

DLtInv = @(sol)(-fluid.dkr(sol)*(1./fluid.mu)'...
                ./((1-sum(fluid.sr))*fluid.Lt(sol).^2));
numW   = numel(W);
numCF  = size(S.B, 1);
obj    = objectiveFunction(param, G, S, W, rock, fluid, simRes, schedule, controls);

k      = 1;
l_v    = secAdjRes(k+1).resSol.cellFlux;
v      = simRes(k+1).resSol.cellFlux;
dim    = DLtInv(simRes(k+1).resSol);
cellNo = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
S.C    = sparse(1:numel(cellNo), cellNo, 1);
Hv     = S.C'*spdiags( (S.C*dim).*v , 0, numCF, numCF)*S.DK*l_v;
Hv     = Hv + obj.partials(k+1).u2;

for wellNr = 1:numW
    w       = W(wellNr);
    q       = simRes(k+1).wellSol(wellNr).flux; % q_w^{n+1}
    l_q     = secAdjRes(k+1).wellSol(wellNr).flux; % lam_q_w^{n+1}
    Hv      = Hv - w.S.C'*( ( (w.S.C*dim).*q ).*(w.S.DK*l_q) );
end
return
                              