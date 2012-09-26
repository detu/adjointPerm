function h = hash(filename)
%Compute MD5 checksum hash of file.
%
% SYNOPSIS:
%   h = hash(file)
%
% PARAMETERS:
%   file - Name (string) of file on which to compute MD5 checksum hash.
%
% RETURNS:
%   h    - MD5 checksum hash of 'file'.
%
% SEE ALSO:
%   md5sum.

%{
#COPYRIGHT#
%}

% $Id: hash.m 1953 2009-03-31 10:54:12Z bska $

  [err, result]=system(['md5sum ', filename]);
  if err
    error(['Unable to obtain hash for file ', filename]);
  else
    h = strtok(result);
  end
end
