function prm = merge_options(prm, varargin)
%Override default control options.
%
% SYNOPSIS:
%   prm = merge_options(prm, 'pn1', pv1, ...)
%
% PARAMETERS:
%   prm       - Original/default control option structure.  The contents of
%               this structure is problem specific and defined by caller.
%
%   'pn1'/pv1 - 'key'/value pairs overriding default options in `prm'.  A
%               WARNING is issued if a given 'key' is not already present
%               in FIELDNAMES(prm).  Key names are case insensitive.
%
%               The ``message identifier'' of this warning is
%
%                   [<FUNCTIONNAME>, ':Option:Unsupported']
%
%               with <FUNCTIONNAME> being the name of the function calling
%               MERGE_OPTIONS or the string 'BASE' if MERGE_OPTIONS is
%               called directly from the base workspace.
%
% EXAMPLE:
%   prm = struct('foo', 1, 'bar', pi, 'baz', true)
%   prm = merge_options(prm, 'foo', 0, 'bar', rand(10), 'FimFoo', @exp)
%
% RETURNS:
%   prm - Modified parameter structure.
%
% SEE ALSO:
%   FIELDNAMES, WARNING, STRUCT.

%{
#COPYRIGHT#
%}

% $Id: merge_options.m 2622 2009-09-02 11:14:57Z steink $

if nargin > 1,
   if mod(numel(varargin), 2) == 0 && ...
      all(cellfun(@ischar, { varargin{1 : 2 : end} })),
      st = dbstack(1);
      try
         caller = st(1).name;
      catch  %#ok
         caller = 'BASE';
      end
      ofn = fieldnames(prm);
      nfn = { varargin{1 : 2 : end} };
      nfv = { varargin{2 : 2 : end} };

      for n = 1 : numel(nfn),
         ix = strmatch(lower(nfn{n}), lower(ofn), 'exact');
         if ~isempty(ix),
            %prm.(nfn{n}) = nfv{n};
            prm.(ofn{ix}) = nfv{n};
         else
            warning([caller, ':Option:Unsupported'], ...
                    ['Option `', nfn{n}, ''' is not supported']);
         end
      end
   else
      error(msgid('Input:Huh'), ...
            'Huh? Did you remember to unpack VARARGIN?!?');
   end
end
