function [x1x2] = xprod( x1, x2, usr_par)
%
% Version June 6, 2008
% Matthias Heinkenschloss
%


% Euclidean inner product
 x1x2 = x1'*x2;
% x1x2 = sum(x1'*x2)/(size(x1,1) + size(x2,1));

% End of xprod
