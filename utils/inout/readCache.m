function success = readCache(arg)
%Read cached variables associated with given file into caller's workspace.
%
% SYNOPSIS:
%   success = readCache(file)
%
% PARAMETERS:
%   file    - Name (string) of file with which to associate a cache.
%
% RETURNS:
%   success - Whether or not cache reading succeeded.  If successful,
%             the cached variables are loaded directly into the caller's
%             workspace.
%
% SEE ALSO:
%   writeCache.

%{
#COPYRIGHT#
%}

% $Id: readCache.m 1953 2009-03-31 10:54:12Z bska $

success = false;
if iscell(arg)

  id = md5sum(arg);
  fn = ['.cache', filesep, id, '.mat'];
  if exist(fn, 'file') == 2,
    dispif(true, ['Reading ', fn, '\n']);
    evalin('caller', ['load ', fn]);
    success = true;
  end



else

  if ~ischar(arg),
    error('readCache takes a string or cell array ar argument');
  end


  [path, name, type] = fileparts(arg);
  fp = fopen(['.cache', filesep, name, type]);
  if fp > 0,
    % Retrieve old hash of arg from file .cache/fn
    oldhash = fgets(fp);
    fclose(fp);

    % Compare old hash to current hash
    if strcmp(hash(arg), oldhash)
      dispif(true, 'Reading .cache%s%s.mat\n', filesep, oldhash);
      evalin('caller', ['load .cache', filesep, oldhash]);
      success = true;
    end
  end
end
