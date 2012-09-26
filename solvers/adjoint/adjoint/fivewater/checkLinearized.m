function [linSimResNum, linSimRes] = checkLinearized()
% check linearized code by comparing to a perturbed control-vector
% u + epsilon*v  giving a perturbed state vector x + epsilon*w
% For small epsilon, w and v should then satify the linearized 
% equations F_x*w = -F_u*v

epsilon = 1e-6;
initSimpleModel

numW = numel(W);
u    = [controls.well.values]';
u    = u(:);
dimU = length(u);

v    = rand(dimU, 1) - .5;

%-------------- Initial and perurbed sim ----------------------------------
simRes  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

u_p     = u + epsilon*v;
controls  = updateControls(controls, u_p);
schedule  = updateSchedule(controls, schedule);
simRes_p  = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);

linSimResNum = simRes;
for k = 2 : controls.numControlSteps+1
    linSimResNum(k).resSol.cellPressure = ...
        (simRes_p(k).resSol.cellPressure-simRes(k).resSol.cellPressure)/epsilon;
    linSimResNum(k).resSol.facePressure = ...
        (simRes_p(k).resSol.facePressure-simRes(k).resSol.facePressure)/epsilon;
    linSimResNum(k).resSol.cellFlux = ...
        (simRes_p(k).resSol.cellFlux-simRes(k).resSol.cellFlux)/epsilon;
    linSimResNum(k).resSol.faceFlux = ...
        (simRes_p(k).resSol.faceFlux-simRes(k).resSol.faceFlux)/epsilon;
    linSimResNum(k).resSol.sw = ...
        (simRes_p(k).resSol.sw-simRes(k).resSol.sw)/epsilon;
    for k1 = 1:numW
        linSimResNum(k).wellSol(k1).flux = ...
            (simRes_p(k).wellSol(k1).flux-simRes(k).wellSol(k1).flux)/epsilon;
        linSimResNum(k).wellSol(k1).pressure = ...
            (simRes_p(k).wellSol(k1).pressure-simRes(k).wellSol(k1).pressure)/epsilon;
    end
end
    

%-------------- Linearized model -----------------------------------------
controls  = updateControls(controls, v);
schedule  = setToZero( schedule );
schedule  = updateSchedule(controls, schedule);
linSimRes = runLinearized(simRes, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);


%-------------- Check output if nargout == 0 ---------------------------

if nargout == 0
    disp(['Relative differences (epsilon = ' num2str(epsilon) ' ):'])
    for k = 2 : controls.numControlSteps+1
        ep = norm(linSimResNum(k).resSol.cellPressure - linSimRes(k).resSol.cellPressure)./...
            norm(linSimResNum(k).resSol.cellPressure);
        ef = norm(linSimResNum(k).resSol.cellFlux - linSimRes(k).resSol.cellFlux)./...
            norm(linSimResNum(k).resSol.cellFlux);
        es = norm(linSimResNum(k).resSol.sw - linSimRes(k).resSol.sw)./...
            norm(linSimResNum(k).resSol.sw);
        disp(['Step ' num2str(k-1) ' : ' num2str(ep) '  ' num2str(ef) '  ' num2str(es)])
    end
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
    
    