function wellSol = packageWellSol(flux, press, fW, hW)
%Convert well fluxes and pressures to structure form.
%
% SYNOPSIS:
%   wellSol = packageWellSol(flux, press, fW, hW)
%
% PARAMETERS:
%   flux   - Vector of well fluxes.  One scalar value for each well
%            *perforation* in the model.  Assumed to be stored as
%
%              flux = [wf_11, wf_12, ..., wf_1n1, ...
%                      wf_21, wf_22, ..., wf_2n2, ...
%                      ...
%                      wf_k1, wf_k2, ..., wf_knk]
%
%            with 'k' being the number of wells, 'wf_ij' being the well
%            flux in perforation 'j' of well 'i', and 'ni' being the total
%            number of perforations in well 'i'.
%
%   press  - Vector of well pressures.  One scalar (non-negative) value for
%            each well in the model.
%
%   fW, hW - Output arrays 'fW' and 'hW' respectively from function
%            'unpackWellSystemComponents'.
%
% RETURNS:
%   wellSol - One-dimensional well solution structure array, one element
%             for each well in the model.  Each element, 'wellSol(i)', is a
%             structure having the following fields:
%                - flux     -- All perforation fluxes for well 'i'.
%                - pressure -- The well pressure in well 'i'.
%
% SEE ALSO:
%   unpackWellSystemComponents.

%{
#COPYRIGHT#
%}

% $Id: packageWellSol.m 1953 2009-03-31 10:54:12Z bska $

   nW = sum(cellfun(@(x) numel(x) > 0, fW));         % Number of wells.

   ix_f = cumsum([0, cellfun(@numel, fW)]);
   ix_p = cumsum([0, cellfun(@numel, hW)]);

   wellSol = repmat(struct('flux', [], 'pressure', []), [1, nW]);
   for w = 1 : nW,
      wellSol(w).flux     = -flux (ix_f(w) + 1 : ix_f(w + 1));
      wellSol(w).pressure =  press(ix_p(w) + 1 : ix_p(w + 1));
   end
end
