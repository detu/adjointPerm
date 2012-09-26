function writeCache(arg, varargin)
% Write callers workspace to matfile in directory ./.cache/
%
% SYNOPSIS:
%   writeCache(fn)
%   writeCache(fn, {'var1', 'var2', ...})
%   writeCache({arg1, ...}, {'var1', 'var2', ...})
%
% PARAMETERS:
%   All variables in callers workspace.
%
% RETURNS:
%   none
%
% SEE ALSO:
%   readCache.

%
% save -append
% save struct field1 field2 ...
%

%{
#COPYRIGHT#
%}

% $Id: writeCache.m 1953 2009-03-31 10:54:12Z bska $

 if nargin == 2
   var = varargin{1};
 else
   var = {};
 end

 if iscell(arg)
   if ~isdir('.cache')
     mkdir('.cache');
   end
   id = md5sum(arg);
   fn = ['.cache', filesep, id];

   dispif(true, 'Saving to %s\n', fn);
   evalin('caller', ['save ', fn, sprintf(' %s', var{:})]);

 else

 [path, name, type]=fileparts(arg);


 % old hash exist and is equal to current hash
 fp = fopen(['.cache/', name, type], 'r');
 h  = hash(arg);

 if ~isdir('.cache')
   mkdir('.cache');

 elseif fp ~= -1
   oldhash = fgets(fp);

   %if hash has not changed
   if strcmp(oldhash, h)
     dispif(true, 'Not writing cache file.');
     return

   %if hash has changed
   else
     dispif(true, 'Removing   .cache/%s.mat\n', oldhash);
     system(['rm .cache/', oldhash,'.mat']);
     system(['rm .cache/', name, type]);
   end
 end

 fp = fopen(['.cache/', name, type], 'w');
 fprintf(fp, '%s', h);
 fclose(fp);

 % Save
 dispif(true, 'Saving to .cache/%s.mat\n', h);
 evalin('caller', ['save .cache/', h, sprintf(' %s',var{:})]);
end
