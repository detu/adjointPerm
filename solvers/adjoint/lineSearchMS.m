function [simRes, schedule, controls, data] = lineSearchMS(...
    simRes, G, CG, SG, S, CS, p, W, rock, fluid, schedule, controls, grad, objectiveFunction, varargin)
% Run ad-hoc line search based on given gradient 
%
% SYNOPSIS:
%  [simRes, schedule, controls, data] = lineSearch(...
%    simRes, G, S, W, rock, fluid, schedule, controls, grad, objectiveFunction, varargin)
%    
% DISCRIPTION:
%  1) The algorithm first checks if any of the elems of u lies on the
%     boundary of the ineq. constr. with elems of grad pointing
%     outwards. If this is the case, the corresponding columns of the 
%     linEqCons are set to zero.
%  2) The gradient is modified by a non-negative projection
%     (I-A*inv(A'*A)A') where the diag elems of I are set to zero for 
%     boundary controls 
%  3) The algorithm checks if u + stepSize*grad is a feasible point 
%     given the ineq. constr. If not, stepSize is scaled such that 
%     u + stepSize*grad lies exactly on the boundary. 
%  4) The objective function value is computed for u + stepSize*grad, and 
%     if this value is an improvement, the algorithm is happy and returns, if 
%     not, a line search is performed until an improvement is found
%
% PARAMETERS:
%   G, S, W, rock, fluid - usual structures
%   simRes               - ...
%   schedule, controls   - ...
%   grad                 - gradien as given by computeGradient
%   objectiveFunction    - handle to objective function
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%          - MaxPoints           : maximum number of points before to
%                                  terminate line search algorithm
%          - StepSize            : Step size 
%          - BoundaryTol         : tollerance at which controls are considered to 
%                                  lay excactly on the boundary of some ineq. constr.
%          - VerboseLevel        : amount of output to screen
%           
% RETURNS:
%   simRes     - Results for 'best' forward simulation
%   schedule   -  
%   controls   - 
%   data       - structure with fields
%                value         : objective function value for best run
%                relNormGrad   : relative norm of gradient |du|/|u|
%                success       : whether or line search procedure succeeded
%                fraction      : optimal fraction obtained during the line
%                                search
%
%
% SEE ALSO:
%   
% 
opt     = struct('MaxPoints',  30, ...
                 'StepSize',        1, ...
                 'BoundaryTol',     0, ...
                 'VerboseLevel',    0);
opt          = merge_options(opt, varargin{:});
numPnts      = opt.MaxPoints;
stepSize     = opt.StepSize;
relBTol      = opt.BoundaryTol;
verboseLevel = opt.VerboseLevel;

%--------------------------------------------------------------------------
if verboseLevel >= 0, fprintf('\n********* Starting line search ***********\n'); end

numSteps        = numel(grad);
numControlWells = numel(controls.well);

uInit    = [controls.well.values]';
uInit    = uInit(:);
normU    = norm(uInit, 'inf');  

obj         = objectiveFunction(G, S, W, rock, fluid, simRes);

%--------------------------------------------------------------------------
% ----------- Scale gradient wrt control-types  ---------------------------
[grad, scaleFacs] = scaleGradient(grad, simRes, controls, obj);

%--------------------------------------------------------------------------
minMax   =  repmat( vertcat(controls.well.minMax), numSteps, 1);

% Find which controls lie on boundary
bControlsMax = ( abs( uInit - minMax(:,2) ) < relBTol * normU );
bControlsMin = ( abs( uInit - minMax(:,1) ) < relBTol * normU );
%bControls    = bControlsMax || bControlsMin;

%uInit(bControlsMax) = maxValues(bControlsMax);
%uInit(bControlsMin) = minValues(bControlsMin);

%--------------------------------------------------------------------------
% -----------Find feasible initial maxStepSize  ---------------------------

du      = cell2mat( {grad{:}}');
gradPos = (du > 0);
gradNeg = (du < 0);

% Set gradient to zero for boundary controls with outward pointing
% gradient => fix controls, project the gradient accordingly
fixedControls = or( and(bControlsMax, gradPos ), ...
                    and(bControlsMin, gradNeg ) );
newFixedControls = fixedControls;
while any( newFixedControls )
    FC = reshape(fixedControls, numControlWells, numSteps);
    for step = 1 : numSteps
        fc      = FC(:, step);
        if any(fc)
            g          = grad{step};
            grad{step} = projection(controls, fc) * g;
        end
    end
    %Recompute
    du      = cell2mat( {grad{:}}');
    gradPos = (du > 0);
    gradNeg = (du < 0);
    newFixedControls = or( and(bControlsMax, gradPos ), ...
                           and(bControlsMin, gradNeg ) );
    fixedControls    = or( fixedControls, newFixedControls);
end

maxStepSizes          = ones( size(uInit) );
maxStepSizes(gradPos) = (minMax(gradPos, 2) - uInit(gradPos))./du(gradPos);
maxStepSizes(gradNeg) = (minMax(gradNeg, 1) - uInit(gradNeg))./du(gradNeg);
maxStepSize           = min(maxStepSizes);

% Set stepSize to maximum allowed
maxStepSize = min(maxStepSize, stepSize);

%--------------------------------------------------------------------------
%------------ Line search: find optimal alph within [0 stepSize] ----------

objValues   = -inf*ones(1, 4);
objValues(1)= obj.val; 


gri         = 2/(1 + sqrt(5));  % inverse golden ration
% points along line x1 --- x2 -- x3 --- x4
x = [0, 1-gri, gri 1]*maxStepSize;

% Compute values for x2,x3,x4
for k = 2 : numPnts
    if k <= 4
        ll   = [0 4 3 2];
        indx = ll(k);
    elseif k > 4
        [indx, x, objValues] = rejectPoint(x, objValues);
    end
    if verboseLevel >= 0, fprintf('\nPoint %2d. Forward simulation: ', k);tic; end
    alph      = x(indx);
    uCur      = uInit + alph*du;
    controls  = updateControls(controls, uCur);
    schedule  = updateSchedule(controls, schedule);
    simRes    = runScheduleMS(simRes(1).resSol, G, CG, SG, S, CS, p, W, rock, fluid, schedule, ...
                            'VerboseLevel', verboseLevel);
    obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
    objValues(indx) = obj.val;
    if verboseLevel >= 0
        tt = toc;
        fprintf('%4.2f sec. Objective value = %6.4f', tt, obj.val); 
    end
    [val, pos] = max(objValues);
    if or(and(k==3, pos == 4), and(k > 9, x(1)~=0) ), break; end
end

%----------------------------------------------------------------------------------------
% Find largest value, and recompute resSim (not neccessary if saved,
% should be done ...)
[val, pos] = max( objValues);
alph       = x(pos);
if alph == 0
    success = false;
    if verboseLevel >= 0, fprintf('\nLine search failed\n');tic; end
else
    if verboseLevel >= 0, fprintf('\nOptimal fraction: %5.4f\n', alph/maxStepSize);tic; end
    success   = true;
    uCur      = uInit + alph*du;
    controls  = updateControls(controls, uCur);
    schedule  = updateSchedule(controls, schedule);
    simRes    = runScheduleMS(simRes(1).resSol, G, CG, SG, S, CS, p, W, rock, fluid, schedule, ...
                            'VerboseLevel', verboseLevel);
end

% Update data-struct
data.value       = val;
data.relNormGrad = norm(du ./ scaleFacs, Inf);
data.success     = success;
data.fraction    = alph/maxStepSize;
data.stepSize    = alph;



function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end

function projector = projection(controls, fixedControls)
I = eye( numel(controls.well) );
I(:, fixedControls) = 0;
projector = I;
if ~isempty( controls.linEqConst )
    A = controls.linEqConst.A;
    A(:, fixedControls) = 0;
    projector = projector - A'*( (A*A')\A );
end

function [indx, x, objValues] = rejectPoint(x, objValues)
gri         = 2/(1 + sqrt(5));
[v, pos] = min( objValues([1 4]) );
if pos == 1         % reject x1
    x           = [x(2) x(3) 0 x(4)];
    x(3)        = x(1) + gri*( x(4) - x(1) );
    objValues   = [objValues(2) objValues(3) 0 objValues(4)];
    indx        = 3;
else
    x           = [x(1) 0 x(2) x(3)];
    x(2)        = x(1) + (1-gri)*( x(4) - x(1) );
    objValues   = [objValues(1) 0 objValues(2) objValues(3)];
    indx        = 2;
end
    
function [grad, scaleFacs] = scaleGradient(grad, simRes, controls, obj)
nSteps = numel(simRes) - 1;
cq = zeros(nSteps, 1); cp = zeros(nSteps, 1);
for step = 1 : nSteps
    ws  = simRes(step+1).wellSol;
    q   = cellfun(@sum, {ws.flux} );
    cq(step) = .5*(max(q) - min(q));
    p   = [ws.pressure];
    cp  = max(p) - min(p);
end
cq = mean(cq); cp = mean(cp);

controlTypes = {controls.well.type}';
isRate       = strcmp( controlTypes , 'rate');
fac          = ( isRate*cq + (~isRate)*cp );

nCSteps = controls.numControlSteps;
for cStep = 1 : nCSteps
    grad{cStep} = ( (fac.^2)/obj.val ) .* grad{cStep};
end

scaleFacs = repmat( fac, nCSteps, 1 );






