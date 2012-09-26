function [simRes] = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, varargin)
% runSchedule -- Run simulation based on schedule.
%
% SYNOPSIS:
%   simRes = runSchedule(resSolInit, G, S, W, fluid, schedule, pn, pv, ...)
%
% DESCRIPTION:
%   
% PARAMETERS:
%   resSolInit  -
%   G           - Grid data structure.
%   S           -
%   W           -
%   fluid       -
%   schedule    -
%
%
% RETURNS:
%   simRes      - (numSteps+1) x 1 structure having fields
%                   - timeInterval
%                   - resSol
%                   - wellSol
%
%
% SEE ALSO:
%  
opt     = struct('Verbose',  false , ...
                 'VerboseLevel', 2);
opt     = merge_options(opt, varargin{:});

verboseLevel2 = opt.Verbose || (opt.VerboseLevel == 2);
verboseLevel1 = opt.Verbose || (opt.VerboseLevel > 0);

numSteps = numel(schedule);
resSol   = resSolInit;
wellSol  = initWellSol(W, 0); 

% Initial conditions
simRes(1).timeInterval  = [0 0];
simRes(1).resSol        = resSol;
simRes(1).wellSol       = [];
if verboseLevel2, dispSchedule(schedule); end

% Pick solver
if strcmp(S.type, 'hybrid')
   solver = 'hybrid';
else
   solver = 'mixed';
end

if verboseLevel1, fprintf('\n******* Starting forward simulation *******\n'); end
for k = 1 : numSteps
    if verboseLevel1, fprintf('Time step %3d of %3d,   ', k, numSteps); end
    W     = updateWells(W, schedule(k));
    interval = schedule(k).timeInterval;
    dt       = interval(2) - interval(1);
    
    % ---- Pressure Equation -----
    if verboseLevel1, fprintf('Pressure:'); tic; end
    [resSol, wellSol] = solveIncompFlow(resSol, wellSol, G, S, fluid, ...
                                        'wells', W, 'Solver', solver);
    
    if verboseLevel1, t = toc; fprintf('%9.3f sec,   ', t); end
    
    % ---- Saturation Equation ---
    if verboseLevel1, fprintf('Transport:'); tic; end
                         
    resSol = implicitTransport(resSol, wellSol, G, dt, rock, fluid,  ... 
                               'wells', W,                           ...
                               'nltol', 1.0e-6, 'lstrials', 50,      ...
                               'maxnewt', 100,  'tsref',  15,        ...
                               'verbose', opt.Verbose);      
                            
    if verboseLevel1, t = toc; fprintf('%9.3f sec\n', t); end                         
    
    % update simRes structure
    simRes(k+1).timeInterval  = interval;
    simRes(k+1).resSol        = resSol;
    simRes(k+1).wellSol       = wellSol;
end

    