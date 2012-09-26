function A = rldecode(A, n, dim)
%Decompress run length encoding of array A along dimension dim.
%
% SYNOPSIS:
%   B = rldecode(A, n, dim)
%   B = rldecode(A, n)       dim assumed to be 1
%
% PARAMETERS:
%   A         - encoded array
%   n         - repetition of each layer along dimension dim
%   dim       - dimension of A where run length encoding is done.
%               dim > 0, dim < ndims(A).
%
% RETURNS:
%   B         - uncompressed A
%
% EXAMPLE:
%   % 1) Numerical example:
%   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
%   [B,n] = rlencode(A,1)
%   C = rldecode(B,n,1)
%
%   % 2) Retrive 'first' column of G.cellFaces (see grid_structure):
%   G = cartGrid([10, 10, 2]);  
%   cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
%   disp(['CellFace nr. 10 belongs to cell nr: ', num2str(cellNo(10))]);
%
% SEE ALSO:
%   rlencode

%{
#COPYRIGHT#
%}

% $Id: rldecode.m 2338 2009-06-05 17:19:30Z bska $

if nargin < 3
  dim = 1;
end

assert ( all( n>=0 ) );
assert ( numel(n) == size(A, dim) );


% Take dimension we compress along to be first dimension,
% i.e., swap dimensions 1 and dim.
d      = 1:ndims(A);
d([1, dim])   = [dim, 1];
B      = permute(A,d);

r      = n~=0;
B      = reshape(B(r, :), sum(r), []);



% Insert repeated layers and permute back
i      = cumsum([1; double(reshape(n(r), [], 1))]);
j      = zeros(i(end)-1,1);
j(i(1:end-1)) = 1;

A      = permute(reshape(B(cumsum(j),:), sum(n), []), d);
