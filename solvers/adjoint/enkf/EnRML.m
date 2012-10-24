           function U = EnRML(D, Y, R, A, A0, beta)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   function U = EnRML(D, Y, R, A, A0, beta)
%
% DESCRIPTION:
%   (global) EnRML update
%
% PARAMETERS:
%   Y           -   ensemble of perturbed observations
%   R           -   array with measurement error variances
%   D           -   ensemble of simulated observations
%   A           -   ensemble of simulated states at current iteration
%   A0          -   ensemble of simulated states at iteration 1
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

           % EnRML - currently only diagonal R
           S = D - repmat(mean(D,2), 1, nmembers);                
           M = A - repmat(mean(A,2), 1, nmembers);
           G = svdbksb(M', S');
           clear M        
           M = (A0 - repmat(mean(A0,2), 1, nmembers)) / sqrt(nmembers-1);
           GM = G' * M;
           if size(GM,1) <= size(GM,2)
             W2 = (GM*GM' + diag(R));
             W  = GM' / W2;
             clear W2
           else
             W = (eye(nmembers) + (GM' / diag(R)) * GM)  \ GM' / diag(R);
           end
           innovation = D - Y - G' * (A - A0);
           MW = M * W;
           % global update
           U = A0 - MW * innovation;
           if beta ~= 1
             U = beta * U + (1-beta) * A; 
           end
           clear W G GM 
