function [x1x2] = xprod1( x1, x2, usr_par)
%
% Version June 6, 2008
% Matthias Heinkenschloss
%


% Euclidean inner product
% x1x2 = x1'*x2;
% x1x2 = sum(x1'*x2)/size(x1,1);
x1x2 = sum(x1.*x2/numel(x1));
% x1x2 = sum(x1.*x2/1e10);

% End of xprod
