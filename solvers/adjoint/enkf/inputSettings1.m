%----------------------------------------------------------------------------------
% DESCRIPTION:
%   Example input file for mrstEnKF. This example: 2D reservoir, 5-spot
%   well pattern. Random Gaussian field are provided to generate ensemble
%   realizations.
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%----------------------------------------------------------------------------------

%----------------------------------------------------------------------------------
% General settings
%----------------------------------------------------------------------------------
 codedir        = ROOTDIR;
 workdir        = [ROOTDIR 'TEST1_HM/'];
 modelname      = 'testGrid1';

%---------------------------------------------------------------------------------- 
% EnKF settings
%----------------------------------------------------------------------------------
 method         = 'EnKF';           % updating method (EnKF, LEnKF, EnRML)
 nmembers       = 50;               % ensemble size
 restart        = 1;                % 0: continue simulation after each update
                                    % 1: restart from start time after each update
 inflation      = 0;                % 0: no inflation, 1: apply after each update
 transform      = {1 0 ...          % transformation functions for the states
                   0 [2 0.2 0.8]};  % (see applyTransform.m)
 niterations    = 1;                % max number of iterations to same update time
                                    % 0: no update, 1: one update, 2: use same
                                    % data twice, etc.
 nrepeats       = 0;                % rumber of repeats of the experiment
 
%---------------------------------------------------------------------------------- 
% Time steps
%----------------------------------------------------------------------------------
 tsim           =  30.00*year();    % simulation end time
 tend           =  20.00*year();    % history match end time
 dt             =   0.50*year();    % simulation time step
 dto            =   2.50*year();    % observation time interval (multiple of dt)
 dtu            =   5.00*year();    % update time interval

%---------------------------------------------------------------------------------- 
% Observations
%----------------------------------------------------------------------------------
 obstype        = {'flx' 'bhp' 'wct'};    % observation type (flx, bhp, wct, sat)
 transformobs   = {0 0 0};                % transformation parameters for observations
 sigma          = [2.5 * meter^3/day, ... % observation error standard deviations
                    25 * barsa(),     ... % for each type
                   0.1 ];
 range          = 2;                      % localization range: c(range) = 0, only
                                          % used with LEnKF
               
%---------------------------------------------------------------------------------- 
% End of input
%----------------------------------------------------------------------------------

 itransform = transform;
 for i = 1 : length(transform)
     itransform{i}(1) = -1 * transform{i}(1);
 end

 obsnum = zeros(1,numel(obstype)); % 1: flx,2: bhp,3: wct,4: sat
 for i = 1:numel(obstype)
     obsnum(i) = find(strncmp(obstype{i},{'flx' 'bhp' 'wct' 'sat'},3),1);
 end
 
 G = []; % initialisation of geo-model parameters
 
 settings = struct(...
     'codedir',     codedir,    ...
     'workdir',     workdir,    ...
     'modelname',   modelname,  ...
     'method',      method,     ...
     'nmembers',    nmembers,   ...
     'restart',     restart,    ...
     'inflation',   inflation,  ...
     'transform',   transform,  ...
     'niterations', niterations,...
     'nrepeats',    nrepeats,   ...
     'tend',        tend,       ...
     'dt',          dt,         ...
     'dto',         dto,        ...
     'dtu',         dtu,        ...
     'obsnum',      obsnum,    ...
     'sigma',       sigma);
 
 if niterations == 0
     restart = 0;
 end
 