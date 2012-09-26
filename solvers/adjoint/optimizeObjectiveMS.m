function [simRes, schedule, controls, output] = ...
                optimizeObjectiveMS(G, CG, SG, S, CS, p, W, rock, fluid, resSolInit, ...
                schedule, controls, objectiveFunction, varargin)
% optimizeObjectiveMS -- Run whole optimization proccess multiscale
%
% SYNOPSIS:
%  [simRes, schedule, controls, output] = optimizeObjectiveMS( ...
%                       G, CG, S, CS, W, rock, fluid, resSolInit, ...
%                       schedule, controls, objectiveFunction, varargin) 
%
% PARAMETERS:
%   G, S, W, rock, fluid - usual structures
%   resSolInit           - initial 'solution' minimum containing field
%                          resSol.sw
%   schedule, controls   - ...
%   objectiveFunction    - handle to objective function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%          - Tol           : tollerance in line search algorithm, proccess stops 
%                            when |du|/|objective| < Tol (default 1e-2)
%          - VerboseLevel  : level of output to screen during progress (default: 0).
%                              < 0  : no output
%                                0  : minumal
%                                1  : quite a bit
%                                2  : loads
%          - PlotProgress  : whether or not to plot output during progress
%                            (default: true)
%          - OutputLevel   : level of output given in structure output (default: 0)
%                              < 0  = nothing (empty)
%                                0  = objective function value for every
%                                     iteration
%                                1  = ?
%           
% RETURNS:
%   simRes     - Results for last (optimal) simulation
%   schedule   -  
%   controls   - 
%   outPut     - as described above
%
% SEE ALSO:
%   

opt     = struct('Tol',           1e-4, ...
                 'VerboseLevel',     0, ...
                 'PlotProgress',  true, ...
                 'OutputLevel',      1);
opt     = merge_options(opt, varargin{:});
tol          = opt.Tol;            % Relative improvement tolerance in objective function
verboseLevel = opt.VerboseLevel;
plotting     = opt.PlotProgress;
outputLevel  = opt.OutputLevel;
stepConst    = 1; % Hack
%--------------------------------------------------------------------------
objValues           = [];
iterationNum        = 0;
relNormGrad         = inf;
successLineSearch   = true;
stepSize            = 1; %stepConst*.1/schedule(1).timeInterval(2);            % Initial stepsize used in line search

schedules{1} = schedule; 

while ( relNormGrad >= tol) && successLineSearch && iterationNum < 30
    iterationNum = iterationNum + 1;
    if verboseLevel >= 0
        fprintf('\n********** STARTING ITERATION %3d ****************\n', iterationNum);
    end
    
    % FORWARD SOLVE
    if iterationNum == 1   % For >1, forward sim is done in lineSearch
        if verboseLevel == 0, fprintf('\nForward solve %3d: ', iterationNum); tic; end;
        
        simRes = runScheduleMS(resSolInit, G, CG, SG, S, CS, p, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);
        
        if verboseLevel == 0, tt = toc; fprintf('%6.3f sec. \n', tt);end
        
        obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
        objValue  = obj.val;
        objValues = objValue;
        
        if verboseLevel >= 0, fprintf('Initial function value: %6.3f\n', objValue); end
        if plotting, plotProgress(G, simRes, controls, schedule, objValues); end
    end


    % ADJOINT SOLVE
    if verboseLevel == 0, fprintf('Adjoint solve %3d: ', iterationNum); tic; end;
    
       adjRes = runAdjointMS(simRes, G, CG, SG, S, CS, W, rock, fluid, schedule, controls, ...
                        objectiveFunction, 'VerboseLevel', verboseLevel);
    
    if verboseLevel == 0, tt = toc; fprintf('%6.3f sec. \n', tt);end
    
    
    % COMPUTE GRADIENT 
    grad   = computeGradient(W, adjRes, schedule, controls);  

    
    % PERFORM LINE SEARCH
    [simRes, schedule, controls, lsData] = lineSearchMS( ...
                     simRes, G, CG, SG, S, CS, p, W, rock, fluid, schedule, controls, grad, ...
                     objectiveFunction, 'VerboseLevel', verboseLevel, 'StepSize', stepSize);
                 
    objValues         = [objValues; lsData.value];
    relNormGrad       = lsData.relNormGrad;
    successLineSearch = lsData.success;
    if lsData.fraction == 1
        stepSize = 2*stepSize;
    else
        stepSize          = lsData.stepSize;
    end
    schedules{iterationNum} = schedule;
  
    if verboseLevel >= 0
            fprintf('\nObtained function value: %6.3f, relative gradient norm: %6.5f\n', objValues(end), relNormGrad);
            fprintf('\nNext stepsize: %6.5f\n', stepSize);
    end
    
    if plotting     
        plotProgress(G, simRes, controls, schedule, objValues); 
    end
end
output.values = objValues;
output.schedules = schedules;

function [] = plotProgress(G, simRes, controls, schedule, val)
setAxis = false;
%maxY    = max([controls.well.maxVal]);
%minY    = min([controls.well.minVal]);
%axLimsY = [minY maxY];
%if any( isinf(axLimsY) ), setAxis = false; end

times   = [schedule.timeInterval]';
for k = 1 : numel(controls.well)
    v          = [controls.well(k).values]';
    V          = [v(1 : end); v(1 : end)];
    data(:, k) = V(:);
    legStr{k}  = ['u_{', num2str(k), '}'];     
end
clf
subplot(1,3,1)
plot(val);

subplot(1,3,2)
plot(times, data);
title(['Objective function value: ', num2str(val(end))]);
legend(legStr)
if setAxis; axis([ [min(times) max(times)] axLimsY]); end

subplot(1,3,3)
plotCellData(G, 1-simRes(end).resSol.s);
caxis([0 1]);axis tight;
drawnow
    


    


    
    
    
    
    
                 
                 
                 