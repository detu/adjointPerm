function [A] = applyTransform(A, nstat, transform)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   state_transformed = applyTransform (state, nstat, transform)
%
% DESCRIPTION:
%   Compute transformation of state variables to improve ensemble distribution
%   charactersistics. Performs the inverse transformations for negative
%   transformation parameter values.
%
% PARAMETERS:
%   state       -   ensemble matrix of state vectors
%   nstat       -   vector with lengths of each state variable
%   transform   -   vector with transformation parameters
%
% RETURNS:
%   Ensemble matrix with transformed state vectors
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 

for i = 1 : length(transform)
  tf = cell2mat(transform(i));
  if tf(1) ~= 0
      var = A(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:);
      
      switch tf(1)
          
        % log transform
        case 1
          var = convertTo(var,milli*darcy);  
          var = log(var);
        case -1
          var = exp(var);
          var = convertFrom(var,milli*darcy);
          
        % logarithmic transform for bounded variables
        case 2
          var(var == tf(2)) = tf(2) + eps;
          var(var == tf(3)) = tf(3) - eps;
          var = log((var - tf(2))./(tf(3)-var)); 
        case -2
          m = (tf(2)+tf(3))/2;
          d = (tf(3)-tf(2))/2;
          var = m + d * (1-1./exp(var)) ./ (1+1./exp(var));
          var(var < tf(2)+eps) = tf(2);
          var(var > tf(3)-eps) = tf(3);

        % local normal-score transform
        case 3   
          for j = 1 : size(var,1)
            a = var(j,:);  
            b = randn(size(a));
            [dum1, idx1] = sort(a); % dum1 = a(idx1)
            [dum2, idx2] = sort(b);
            for k = 1 : numel(a)
                a(idx1(k)) = dum2(k);
            end
            var(j,:) = a;
            s1(j,:) = dum1;
            s2(j,:) = dum2;
          end
          save('transform.mat','s1','s2');        
        case -3      
          load('transform.mat','s1','s2');
          for j = 1 : size(var,1)
            b = var(j,:);  
            n = numel(b);
            for k = 1 : n
              if b(k) > s2(j,n)
                 % linear extrapolation
                 a(k) = s1(j,n)+(b(k)-s2(j,n))*(s1(j,n)-s1(j,n-1))/(s2(j,n)-s2(j,n-1));
              elseif b(k) < s2(j,1)
                 % linear extrapolation 
                 a(k) = s1(j,1)+(b(k)-s2(j,1))*(s1(j,2)-s1(j,1))/(s2(j,2)-s2(j,1));
              else
                 for l = 1 : n - 1
                    if b(k) >= s2(j,l) && b(k) < s2(j,l+1)
                       % linear interpolation between 2 values
                       a(k) = s1(j,l)+(b(k)-s2(j,l))*(s1(j,l+1)-s1(j,l))/(s2(j,l+1)-s2(j,l));
                    end
                 end
              end
            end         
            var(j,:) = a;
          end
          
        % specify other transformations here
        
      end %switch
      
      A(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:) = var;
      clear var
  end
end