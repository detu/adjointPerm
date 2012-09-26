function schedule = updateSchedule(controls, schedule)
% Update schedule based on controls
numSteps     = numel(schedule);
controlWells = [controls.well.wellNum]';
U = [controls.well.values]';
for k = 1:numSteps
    schedule(k).values(controlWells) = U(:, k);
end


