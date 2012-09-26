function S = setBCSB(S,BC,Dofs,G)
% setBCSB - Chooses which boundary conditions are to be implemented into
%           the system S and calls the correct setBCSB_* function 
%
% SYNOPSIS:
%   S = setBCSB(S,BC,Dofs,G) 
%
% PARAMETERS:
%   S          - The mixed Stokes-Brinkman system
%
%   BC         - Boundary condition structure from a prior call to 'addBCSB' or
%                an empty structure, which will automatically give no-flow conditions.
%
%   Dofs       - The degrees-of-freedom structure
%
%   G          - The grid structure
%
% RETURNS:
%   S   - The updated system structure with boundary conditions implemented

  error(nargchk(4, 4, nargin, 'struct'));
  
%  fprintf('Implementing boundary conditions:\n '); 
 
  % Checks for normal velocity boundary conditions and calls setBCSB_vn
  isVeln = BC.type == 1; %strcmp('velocity_n',BC.type);  
  if find(isVeln)
    % temp: - allow only no-flow velocity 
    assert(all(BC.value(isVeln)==0));
    if ~strcmp(S.element,'TH')    
        S = setBCSB_vn(S,BC);
    end
  end

  % Checks for pressure boundary conditions and calls setBCSB_p
  isPress = BC.type == 2; %strcmp('pressure',BC.type);
  if find(isPress)
    S = setBCSB_p(S,BC,Dofs,G);
  end
  
  % Checks for tangential velocity boundary conditions and calls setBCSB_vt
  isVelt = BC.type == 3; %strcmp('velocity_t',BC.type);
  if find(isVelt)
    S = setBCSB_vt(S,BC);
  end;
  
  
  % Checks for pressure boundary conditions, and if not calls setBCSB_p_unique
  if ~isPress
    BC.p.unique = [opt.regNode,opt.regValue];
    S = setBCSB_p_unique(S,BC);
  end
