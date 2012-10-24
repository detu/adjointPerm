      function [Y, R, D, loc] = makeObservations(wellData,wellDataE,timeSteps,...
                                obsTimes,obstype,transform,sigma,sat,sate)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [Y, R, D, loc] = makeObservations(wellData,wellDataE,timeSteps,...
%                    obsTimes,obstype,transform,sigma) 
%   [Y, R, D, loc] = makeObservations(wellData,wellDataE,timeSteps,...
%                    obsTimes,obstype,transform,sigma,sat,sate)
%
% DESCRIPTION:
%   Generate arrays with actual and simulated observations.
%
% PARAMETERS:
%   wellData    -   true production history
%   wellDataE   -   ensemble simulated production history
%   timeSteps   -   time steps associated with production history
%   obsTimes    -   times from which observations are to be collected
%   obstype     -   measurement type
%   transform   -   desired transformation for each measurement type
%   sigma       -   measurement error standard deviations
%   
%   sat         -   true saturation field
%   sate        -   ensemble of simulated saturation fields
%
% RETURNS:
%   Y           -   ensemble of perturbed observations
%   R           -   array with measurement error variances
%   D           -   ensemble of simulated observations
%   loc         -   array with grid indices for each measurement
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
                            
      nmembers = size(wellDataE(1).flx,2);
      nwells   = length(wellData);
      ntimes   = length(timeSteps);
      ntimeso  = length(obsTimes);
      [~, ~, iobs] = intersect(obsTimes, timeSteps); clear *dummy
      
      Y = []; R = []; D = []; loc = [];
      if ~isempty(iobs)
      
      %------------------------------------------------------------------------------
      % well observations (rates and/or pressures)
      %------------------------------------------------------------------------------
      match=0;
      for i=1:length(obstype)
        if ~isempty(strmatch(obstype{i},{'flx' 'bhp' 'wct'})), match=match+1; end
      end
      if match > 0 && timeSteps(end) > 0
                
        flx = nan*zeros(nwells, ntimes);
        bhp = nan*zeros(nwells, ntimes);
        wct = nan*zeros(nwells, ntimes);
        for i = 1 : nwells
          flx(i,:) = wellData(i).flx;
          bhp(i,:) = wellData(i).bhp;
          wct(i,:) = wellData(i).ffl;
        end
        flx = flx(:,iobs);
        bhp = bhp(:,iobs);
        wct = wct(:,iobs);
        
        flxe = nan*zeros(nwells, ntimes, nmembers);
        bhpe = nan*zeros(nwells, ntimes, nmembers);
        wcte = nan*zeros(nwells, ntimes, nmembers);
        if nmembers > 0
            for i = 1 : nwells
                flxe(i,:,:) = wellDataE(i).flx; % ntimes x nmembers
                bhpe(i,:,:) = wellDataE(i).bhp;
                wcte(i,:,:) = wellDataE(i).ffl;
            end
        end
        flxe = flxe(:,iobs,:);
        bhpe = bhpe(:,iobs,:);
        wcte = wcte(:,iobs,:);
               
        % well locations (currently based on first perforated layer only)
        wloc = nan*ones(size(flx));
        for i = 1 : nwells
            wloc(i,:) = wellData(i).loc(1);
        end
        
        wloc = reshape(wloc, nwells * ntimeso, 1);
        flx  = reshape(flx,  nwells * ntimeso, 1);
        bhp  = reshape(bhp,  nwells * ntimeso, 1);
        wct  = reshape(wct,  nwells * ntimeso, 1);
        flxe = reshape(flxe, nwells * ntimeso, nmembers);
        bhpe = reshape(bhpe, nwells * ntimeso, nmembers);
        wcte = reshape(wcte, nwells * ntimeso, nmembers);
           
      end
      
      %------------------------------------------------------------------------------
      % construct measurement vector
      %------------------------------------------------------------------------------      
%     if iter == 1
        obs = []; sig = []; loc = []; scl = [];
        for i = 1:length(obstype)
          dummy = eval(obstype{i});
          if ~isempty(transform) && length(transform) >= i
            dummy = applyTransform(dummy,size(dummy,1),transform(i));
          end
          scl = [scl; max(dummy)*ones(size(dummy))];
          obs = [obs; dummy];
          sig = [sig; sigma(i)*ones(size(dummy))];
          if ~isempty(strmatch(obstype{i},{'flx' 'bhp' 'wct'}))
              loc = [loc; wloc];
          end
          if ~isempty(strmatch(obstype{i},{'sat'}))
             loc = [loc; -1*(1:length(dummy))'];
          end
          clear dummy
        end
%     end
        
      if nmembers > 0
        obse = [];    
        for i = 1:length(obstype)
          dummy = eval([obstype{i} 'e']);
          if ~isempty(transform) && length(transform) >= i
            dummy = applyTransform(dummy,size(dummy,1),transform(i));
          end
          obse = [obse; dummy]; 
          clear dummy
        end
      end

      clear bhp* flx* wct* sat* wloc
             
      %------------------------------------------------------------------------------                   
      % remove obsolete observations
      %------------------------------------------------------------------------------
      if nmembers > 0
        remove1 = find(max(abs(repmat(obs,1,nmembers) - obse),[],2) < eps);
        remove2 = find(std(obse,[],2) < eps);
        obse(union(remove1,remove2),:) = [];
        obs(union(remove1,remove2)) = [];
        sig(union(remove1,remove2)) = [];
        loc(union(remove1,remove2)) = [];
        scl(union(remove1,remove2)) = [];
      end
      
      %------------------------------------------------------------------------------                   
      % make ensemble of observations by random perturbation
      %------------------------------------------------------------------------------
      if ~isempty(obs) % && iter == 1
        [Y, R] = perturbObservations(obs, sig, nmembers);
      end

      %------------------------------------------------------------------------------                         
      % ensemble of simulated observations
      %------------------------------------------------------------------------------                       
      if nmembers > 0
        D = obse;
      else
        D = [];
      end
      
      % scaling
      Y = Y ./ repmat(scl,1,nmembers);
      D = D ./ repmat(scl,1,nmembers);
      R = R ./ (scl.^2);
      
      clear obs obse 
      
      end