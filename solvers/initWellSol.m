function wellSol = initWellSol(W, p0)
%Initialize well solution data structure.
%
% SYNOPSIS:
%   wellSol = initWellSol(W, p0)
%
% DESCRIPTION:
%   Initialize the well solution structure to uniform well bottom-hole
%   pressures and all-zero well rates.
%
% PARAMETERS:
%   W  - Well data structure as defined by addWell &c.
%   p0 - Initial uniform well bottom-hole pressure (scalar).
%
% RETURNS:
%   wellSol - Initialized reservoir solution structure having fields
%               - flux     -- Well rates in all perforations (== 0).
%               - pressure -- Well bottom-hole pressure (== p0).
%
% SEE ALSO:
%   initResSol, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Id: initWellSol.m 1953 2009-03-31 10:54:12Z bska $

nW = numel(W);
wellSol = repmat(struct('flux', [], 'pressure', []), [1, nW]);

for w = 1 : nW,
   wellSol(w).flux     = zeros([numel(W(w).cells), 1]);
   wellSol(w).pressure = repmat(p0, [numel(W(w).cells), 1]);
end
