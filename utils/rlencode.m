function [A,n] = rlencode(A, dim)
%Compute run length encoding of array A along dimension dim.
%
% SYNOPSIS:
%   [A,n] = rlencode(A, dim)
%
% PARAMETERS:
%   A         - Array
%   dim       - dimension of A where run length encoding is done.
%               dim > 0, dim < ndims(A).
%
% RETURNS:
%   A         - Compressed A where repeated layers are removed.
%   n         - repetition count of repeated layers in original A.
%
% EXAMPLE:
%
%   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
%   [A,n] = rlencode(A,1)
%
%
% SEE ALSO:
%   rldecode

%{
#COPYRIGHT#
%}

% $Id: rlencode.m 1953 2009-03-31 10:54:12Z bska $

if nargin < 2
  dim = 1;
end


% Take dimension we compress along to be first dimension,
% i.e., swap dimensions 1 and dim.

d      = 1:ndims(A);
d(1)   = dim;
d(dim) = 1;
B      = permute(A,d);

% Pick out (1:end,:,...) and (2:end,:,...) in
% a multidimensional way


% Find positions where layers differ
i    = [find(any(B(1:end-1,:) ~= B(2:end,:), 2)); size(B,1)];

% compare differences in position to find run length.
n = diff([0;i]);

% swap dimensions 1 and dim.
sz = size(B);sz(1)=numel(i);
A  = permute(reshape(B(i,:), sz),d);
