function ticif(bool)
%Evaluate function TIC if input is true.
%
% SYNOPSIS:
%   ticif(bool)
%
% PARAMETERS:
%   bool - Boolean variable (true/false).
%
% COMMENTS:
%   Function used for making code cleaner where verbose option is used.
%
% SEE ALSO:
%   tic, tocif, dispif.

%{
#COPYRIGHT#
%}

% $Id: ticif.m 1953 2009-03-31 10:54:12Z bska $

if bool, tic, end
