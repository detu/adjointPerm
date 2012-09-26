function x = mcolon(lo, hi, s)
%Compute vector of consecutive indices from separate lower/upper bounds.
%
% SYNOPSIS:
%   ind = mcolon(lo, hi)
%   ind = mcolon(lo, hi, stride)
%
% PARAMETERS:
%   lo  - Vector of start values (lower bounds).
%   hi  - Vector of end values (upper bounds).
%   s   - Vector of strides.
%         Optional.  Default value: s = ones([numel(lo), 1]) (unit stride).
%
% RETURNS:
%   ind - [lo(1):hi(1)     , lo(2):hi(2)     , ..., lo(end):hi(end)       ]
%   ind - [lo(1):s(1):hi(1), lo(2):s(2):hi(2), ..., lo(end):s(end):hi(end)]
%
% EXAMPLE:
%   lo  = [1 1 1 1]; hi = [2 3 4 5];
%   ind = mcolon(lo, hi)
%
% SEE ALSO:
%   ARRAYFUN.

%{
#COPYRIGHT#
%}

% $Id: mcolon.m 2989 2009-10-13 13:22:29Z bska $

if numel(lo) ~= numel(hi)
  error('In mcolon: lo and hi must have same number of elements!');
elseif numel(lo) == 0,
  x = [];
elseif nargin < 3
   
  % Remove lo-hi pairs where numel(lo:hi)==0
  i    = hi>=lo;
  hi   = hi(i);
  lo   = lo(i);
  
  m    = numel(lo);
  d    = double(hi - lo + 1);
  n    = sum(d);

  x    = ones(1,n);
  x(1) = lo(1);
  x(1+cumsum(d(1:end-1))) = lo(2:m) - hi(1:m-1);
  x    = cumsum(x);
else
   
  % Remove lo-hi-s triplets where numel(lo:s:hi)==0
  i    = hi>=lo & s>0 | hi<=lo & s<0;
  hi   = hi(i);
  lo   = lo(i);
  s    = s(i);
  
  % Compute lo + (0:(hi-lo)/stride)*stride
  hi = floor((hi-lo)./s);


  m    = numel(lo);
  d    = double(hi + 1);
  n    = sum(d);

  % Expand lo  to [lo(1) lo(1) ... lo(2) lo(2) ... lo(end)]
  LO   = zeros(1,n);
  LO(1)=lo(1);
  LO(1+cumsum(d(1:end-1))) = lo(2:m) - lo(1:m-1);
  LO   = cumsum(LO);

  % Expand stride
  S    = zeros(1,n);
  S(1) = s(1);
  S(1+cumsum(d(1:end-1))) = s(2:m) - s(1:m-1);
  S    = cumsum(S);

  x    = ones(1,n);
  x(1) = 0;
  x(1+cumsum(d(1:end-1))) = -hi(1:m-1);


  x    = cumsum(x).*S+LO;
end
