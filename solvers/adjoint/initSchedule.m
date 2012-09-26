function [schedule] = initSchedule(W, varargin)
% initSchedule -- Initialize schedule structure based on well W.
%
% SYNOPSIS:
%   shedule = initSchedule(W, 'pn', pv, ...)
%
% DESCRIPTION:
%   Initialize schedule 
%
% PARAMETERS:
%   W       - well structure
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - NumSteps  :  number of simulation time steps (default 1)
%               - TotalTime :  total simualtion time (default 1)
%               - TimeSteps :  endtime for each time step assuming t_0 = 0 (alterative to two previous pns)
%               - Verbose   :  display schedule with dispSchedule
%
% RETURNS:
%   schedule   - Initialized numSteps x 1 rate schedule structure having fields
%                   - timeInterval     -- [startTime endTime]
%                   - names            -- {name_1, ..., name_n} 
%                   - types            -- {welltype_1, ... , welltype_n}
%                   - values           -- {val_1, ..., val_2}
%                   
%
% SEE ALSO:
%              dispSchedule
%   

opt = struct('Verbose',   false, ...
             'NumSteps',  1, ... 
             'TotalTime', 1, ...
             'TimeSteps', []);
opt = merge_options(opt, varargin{:});

verbose      = opt.Verbose;
timeSteps    = opt.TimeSteps;

if isempty(timeSteps)
    if ~isfield(opt, 'NumSteps') && ~isfield(opt, 'TotalTime')
        error('Non compatible input ...')
    else
        numSteps = opt.NumSteps;
        totTime  = opt.TotalTime;
        timeSteps = (totTime/numSteps)*(1 : numSteps);
    end
end

ts = timeSteps(:);
intervals = [ [0; ts(1:end-1)] ts(1:end) ];

for k = 1 : length(timeSteps)
    schedule(k).timeInterval    = intervals(k, :);
    schedule(k).names           = {W(:).name}';
    schedule(k).types           = {W(:).type}';
    schedule(k).values          = [W(:).val]';
end

if verbose, dispSchedule(schedule); end