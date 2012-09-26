function s = msgid(s)
%Construct Error/Warning message ID by prepending function name.
%
% SYNOPSIS:
%   s = msgid(s)
%
% PARAMETERS:
%   s - A string of the form '[<component>]:<mnemonic>' suitable as an
%       ERROR or WARNING-type message identifier.
%
% RETURNS:
%   s - The same string, though with the name of the function calling
%       'msgid' prepended to 's' (or the string 'BASE' if function 'msgid'
%       is called from the base workspace).
%
% SEE ALSO:
%   error, warning.

%{
#COPYRIGHT#
%}

% $Id: msgid.m 1953 2009-03-31 10:54:12Z bska $

   st = dbstack(1);

   try
      caller = st(1).name;
   catch %#ok
      caller = 'BASE';
   end

   s = [caller, ':', s];
end
