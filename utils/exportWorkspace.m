function exportWorkspace()
%Export caller's workspace to 'BASE' workspace.
%
% SYNOPSIS:
%   exportWorkspace()
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   Nothing, but caller's workspace (all accessible variables) is exported
%   to the 'BASE' MATLAB workspace.
%
% SEE ALSO:
%   assignin, evalin.

%{
#COPYRIGHT#
%}

% $Id: exportWorkspace.m 2773 2009-09-21 08:31:17Z jrn $

    s = evalin('caller', 'whos');
    
    for i = 1 : numel(s),
        cmd = sprintf('assignin(''base'', ''%s'', eval(''%s''))', ...
                      s(i).name, s(i).name);
        try 
           evalin('caller', cmd);
        end
    end
end
