       function [Y, R] = perturbObservations(obs, sig, nmembers)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [Y, R] = perturbObservations(obs, sig, nmembers)
%
% DESCRIPTION:
%   Add an ensemble of perturbations to the observations based on provided 
%   measurement error standard deviations.
%
% PARAMETERS:
%   obs             -   measurements
%   sig             -   measurement error standard deviations 
%   nmembers        -   ensemble size
%
% RETURNS:
%   Y           -   ensemble matrix of perturbed measurement vectors
%   R           -   array with measurement error variances
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%----------------------------------------------------------------------------------
       % observations - currently only diagonal R
       y = obs;
       
       nobs = length(y);
           
       noise = randn(max(10000,nobs),1);
        
       for i = 1 : length(y)
           y(i) = y(i) + sig(i)*noise(end-nobs+i);
       end
       R = sig.^2;
       clear sig noise

       % ensemble of perturbed observations
       if nmembers > 0
           Y = repmat(y, 1, nmembers);
           for i = 1:size(Y,1)
             rndm(i,:) = randn(1,nmembers); 
             rndm(i,:) = rndm(i,:) - mean(rndm(i,:)); 
             rndm(i,:) = rndm(i,:) / std(rndm(i,:));
             Y(i,:) = Y(i,:) + sqrt(R(i)) * rndm(i,:);
           end
           clear rndm
       else
           Y = [];
       end
