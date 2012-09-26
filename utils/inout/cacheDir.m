function tmp = cacheDir
%Define platform-dependent directory for temporary files.

%{
#COPYRIGHT#
%}

% $Id: cacheDir.m 1953 2009-03-31 10:54:12Z bska $

if isunix,
  tmp = '/tmp/';
else
  tmp = ['.', filesep];
end

tmp = [tmp, 'cache', filesep];
