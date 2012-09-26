function nc = dispif(bool, varargin)
%Display value if input is true.
%
% SYNOPSIS:
%   nc = dispif(bool, format, arg, ...);
%
% PARAMETERS:
%   bool   - Boolean variable
%
%   format - SPRINTF format specification.
%
%   arg    - OPTIONAL arguments to complete 'format'.
%
% RETURNS:
%   nc     - Number of characters printed to output device.  If 'bool' is
%            FALSE, then nc=0.
%
% COMMENTS:
%   Function used for making code cleaner where 'verbose' option is used.
%
% SEE ALSO:
%   SPRINTF, tocif.

%{
#COPYRIGHT#
%}

% $Id: dispif.m 2338 2009-06-05 17:19:30Z bska $

nc = 0;
if bool, nc = fprintf(varargin{:}); end
