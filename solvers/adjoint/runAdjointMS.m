function [adjRes] = runAdjointMS(simRes, G, CG, SG, S, CS, W, rock, fluid, schedule, controls, objectiveFunction, varargin)
% runAdjoint -- Run adjoint simulation based on simRes and schedule.
%
% SYNOPSIS:
%   adjRes = runAdjoint(simRes, G, S, W, fluid, schedule, objective, pn, pv, ...)
%
% DESCRIPTION:
%   
% PARAMETERS:
%   simRes      -
%   G           - Grid data structure.
%   S           -
%   W           -
%   fluid       -
%   schedule    -
%   objective   - function handle
%
%
% RETURNS:
%   adjRes      - numSteps x 1 structure having fields
%                   - timeInterval
%                   - resSol
%                   - wellSol
%
%
% SEE ALSO:
%  
opt     = struct('Verbose',  false , ...
                 'VerboseLevel', 0);
opt     = merge_options(opt, varargin{:});
verboseLevel2 = opt.Verbose || (opt.VerboseLevel == 2);
verboseLevel1 = opt.Verbose || (opt.VerboseLevel > 0);

numSteps = numel(schedule);
adjRes   = [];
obj      = objectiveFunction(G, S, W, rock, fluid, simRes, schedule, controls);
if verboseLevel1, fprintf('\n******* Starting adjoint simulation *******\n'); end
for k = numSteps : -1 : 1
    if verboseLevel1, fprintf('Time step %3d of %3d,   ', k, numSteps); end
    W    = updateWells(W, schedule(k));
    
     % ---- Analogue of Saturation Equation 
    if verboseLevel1, fprintf('Transport:'); tic; end
    adjRes = solveAdjointTransportSystemMS(G, SG, S, W, rock, fluid, ...
                                           simRes, adjRes, obj);
    if verboseLevel1, t = toc; fprintf('%9.3f sec,   ', t); end
    
    % ---- Analogue of Pressure Equation
    if verboseLevel1, fprintf('Pressure:'); tic; end
    adjRes = solveAdjointPressureSystemMS(G, CG, S, CS, W, rock, fluid, ...
                                          simRes, adjRes, obj);
    if verboseLevel1, t = toc; fprintf('%9.3f sec\n', t); end
end
