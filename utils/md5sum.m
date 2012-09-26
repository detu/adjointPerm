% md5sum - Compute md5 check sum of all input arguments
%
% SYNOPSIS:
%   str = md5sum(args...)
%
% PARAMETERS:
%   md5sum can take an arbitrary number of arguments.
%
% RETURNS:
%   string of 32 characters with the hexadecimal md5 checksum of
%   all numeric and character arrays that are found as plain arrays,
%   in structs or cell arrays.
%
% EXAMPLE:
%   C{1}=struct('a',1,'b',2);C{3}=speye(4);
%   sum = md5sum(C)
%
% NOTE:
%   This utility is written in C.  It must be compiled with
%
%      mex md5sum.c
%
%   before use.  On older systems, the command is
%
%      med -DOLDMATLAB md5sum.c
%
% SEE ALSO:

% $Id: md5sum.m 1144 2009-01-05 18:04:33Z bska $
