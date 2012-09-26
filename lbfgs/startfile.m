% clear all
% 

% Dimensions
n = 49;
m = 49;

% 0bjective function
J   = [];
nii = 0;

% Gaussian basis vectiors
load structure.mat
[ncells,npp] = size(S);
S = reshape(S,ncells*npp,1);

save structure.dat S -ASCII
save start.mat n m nii npp
save objfun.mat J