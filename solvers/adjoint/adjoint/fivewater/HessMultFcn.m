function w = HessMultFcn(U, lambda, v, auxdata)

[G, S, W, rock, fluid, resSolInit, schedule, controls, objectiveFunction] = deal(auxdata{:});

% % Objective function NPV
% controls  = updateControls(controls, U);
% schedule  = updateSchedule(controls, schedule);
% 
% % FORWARD SOLVE
% simRes    = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
% 
% % ADJOINT SOLVE
% adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls,objectiveFunction);

load simAdj;

% COMPUTE HESSIAN 
% Objective function NPV
controls  = updateControls(controls, v);
schedule  = setToZero(schedule);
schedule  = updateSchedule(controls, schedule);

linSimRes = runLinearized(simRes, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

secAdjRes = runSecondOrderAdjoint(simRes, adjRes, linSimRes, G, S, W, rock, fluid, ...
                                  schedule, controls, objectiveFunction, 'VerboseLevel', 0);

w  = -compHv(W, secAdjRes, schedule, controls);                              
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


