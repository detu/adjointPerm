function [Hv] = Hessvec(dv, u, usr_par)
% function [Hv] = Hessvec(G, S, W, rock, fluid, simRes, adjRes, schedule, controls, objectiveFunction, varargin)

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
% opt     = struct('Verbose',  false , ...
%                  'Vector', []);
% opt     = merge_options(opt, varargin{:});
% if opt.Verbose
%     verboseLevel = 1; 
% else
%     verboseLevel = 0;
% end
% 
[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction] = deal(usr_par{:});
% 
% if ~isempty(opt.Vector)
%     dispif(verboseLevel==1, 'Uptating schedule to incorporate vector'); 
%     controls  = updateControls(controls, opt.Vector);
%     schedule  = setToZero(schedule);
%     schedule  = updateSchedule(controls, schedule);
% end

% FORWARD SOLVE
simRes    = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
% numSteps  = numel(schedule);

% ADJOINT SOLVE
adjRes    = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls,objectiveFunction);

% Objective function NPV
controls  = updateControls(controls, dv);
schedule  = setToZero(schedule);
schedule  = updateSchedule(controls, schedule);

linSimRes = runLinearized(simRes, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

secAdjRes = runSecondOrderAdjoint(simRes, adjRes, linSimRes, G, S, W, rock, fluid, ...
                                  schedule, controls, objectiveFunction, 'VerboseLevel', 0);

Hv  = -compHv(W, secAdjRes, schedule, controls);                              
return
                              
%--------------------------------------------------------------------------

function Hv = compHv(W, secAdjRes, schedule, controls)
bhpWells  = find( strcmp('bhp', {W.type}) );
rateWells = find( strcmp('rate', {W.type}) );
S         = [ W.S ];
Dw        = blkdiag( S.D );
DwD       = Dw(:, bhpWells);

[A_N, b_N, A_D, b_D] = controls2Wells(W, schedule, controls);
Hv = [];
for k = 1 : numel(A_N)
    adjWellPres = [secAdjRes(k+1).wellSol.pressure]';
    l_p         = adjWellPres(rateWells);
    l_q         = vertcat(secAdjRes(k+1).wellSol.flux);
        
    Hv  = [Hv; A_N{k}'*l_p + A_D{k}'*DwD'*l_q];   % non-projected gradient
end
                              
function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end 

function schedule = setToZero(schedule)
% Set all vals in schedule to zero
dims = size( schedule(1).values );
for k = 1 : numel( schedule )
    schedule(k).values = zeros( dims );
end