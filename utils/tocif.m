function [] = tocif(bool)
%Evaluate function TOC if input is true.
%
% SYNOPSIS:
%   tocif(bool)
%
% PARAMETERS:
%   bool - Boolean variable (true/false).
%
% COMMENTS:
%   Function used for making code cleaner where verbose option is used.
%
% SEE ALSO:
%   toc, ticif, dispif.

%{
#COPYRIGHT#
%}

% $Id: tocif.m 1953 2009-03-31 10:54:12Z bska $

if bool, toc, end
