function S = setBCSB_p_unique(S,BC)
% setBCSB_p_unique - Sets the pressure zero level in the system S
%
% SYNOPSIS:
%   S = setBCSB_p_unique(S,BC) 
%
% PARAMETERS:
%   S          - The mixed Stokes-Brinkman system
%
%   BC         - Boundary condition structure from a prior call to 'addBCSB' or
%                an empty structure, which will automatically give no-flow conditions.
%
% RETURNS:
%   S   - The updated system structure with boundary conditions implemented

  error(nargchk(2, 2, nargin, 'struct'));

  if size(BC.p.unique,1) ~= 1,
    error('You can not set a unique pressure in more than one node');
  end
    
  S.P0(BC.p.unique(1),BC.p.unique(1)) = S.P0(BC.p.unique(1),BC.p.unique(1))+1;
  S.Q(BC.p.unique(1))                 = S.Q(BC.p.unique(1))-BC.p.unique(2);
  

