           function U = EnKF (D, Y, R, A, beta)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   function U = EnKF (D, Y, R, A, beta)
%
% DESCRIPTION:
%   (global) EnKF update
%
% PARAMETERS:
%   Y           -   ensemble of perturbed observations
%   R           -   array with measurement error variances
%   D           -   ensemble of simulated observations
%   A           -   ensemble of simulated states
%   beta        -   step size
%
% RETURNS:
%   Updated ensemble of model states U
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
           nmembers = size(A,2);

           % EnKF - currently only diagonal R
           S = (D-repmat(mean(D,2),1,nmembers))/sqrt(nmembers-1);           
           if size(S,1) <= size(D, 2),
             W2 = S*S' + diag(R);
             W = S' / W2;
           else
             W = (eye(nmembers) + (S'/diag(R))*S)  \ S'/diag(R);
           end
           U = A + beta * sqrt(nmembers-1).\(A - repmat(mean(A,2),1,nmembers))*(W*(Y - D));
           
           % diagnostic: updated simulated data
%          Ud  = D + S * (W * (Y - D));  

