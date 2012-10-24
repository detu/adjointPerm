function [X] = svbksb(A,B)
%----------------------------------------------------------------------------------
% SYNOPSIS:
%   [X] = svbksb(A,B))
%
% DESCRIPTION:
% Solves A*X=B for X by SVD decomposition and back substitution. See
% Numerical Recipes.
%
% PARAMETERS:
%   A           -   normal matrix
%   B           -   right-hand side
%   
%
% RETURNS:
%   Y           -   solution
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%----------------------------------------------------------------------------------


m=size(A,1);
n=size(A,2);

[u,y,v]=svd(A);         %A=u*w*v'

for j=1:min(m,n), w(j)=y(j,j); end
w(w<1.e-6*max(w))=0;    %threshold
if m<n, w=[w zeros(1,n-m)]; end
    
for k=1:size(B,2)
for j=1:n
    if w(j)~=0
        tmp(j)=(u(:,j)'*B(:,k))/w(j);
    else
        tmp(j)=0;
    end
end
for j=1:n
    x(j)=v(j,:)*tmp';
end
X(:,k)=x;
end