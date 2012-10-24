%----------------------------------------------------------------------------------
% SYNOPSIS:
%
% DESCRIPTION:
%   Simple inflation routine that removes ensemble mean, multiplies the
%   anomalies to restore the original spread, and adds back the mean.
%
% PARAMETERS:
%
% RETURNS:
%   Inflated ensemble U
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
       for i=1:size(A0,1)
         tmp = A0(i,:); m = mean(tmp); tmp = tmp - m; Ea = std(tmp);
         tmp = U(i,:); m = mean(tmp); tmp = tmp - m; Eu = std(tmp);         
         if Eu > 0 && Ea > 0
           % re-inflate to original spread  
           tmp=tmp*Ea/Eu;
%          % re-inflate only if decrease is significant
%          if sqrt(Eu/Ea) < 0.6 
%            tmp=tmp*Ea/Eu;
%          end
         else
           % add original anomalies  
           tmp=A0(i,:); tmp = tmp - mean(tmp);  
         end
         U(i,:)=tmp+m;
       end
       clear tmp m Eu Ea
