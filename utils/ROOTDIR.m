function root = ROOTDIR
%Retrieve full path of Toolbox installation directory.
%
% SYNOPSIS:
%   root = ROOTDIR
%
% PARAMETERS:
%   none.
%
% RETURNS:
%   root - Full path to MRST installation directory.

%{
#COPYRIGHT#
%}

% $Id: ROOTDIR.m 1953 2009-03-31 10:54:12Z bska $

nm = mfilename('fullpath');
ix = strfind(nm, filesep);
if ~isempty(ix),
   root = nm(1 : ix(end-1));
else
   root = nm;
end
