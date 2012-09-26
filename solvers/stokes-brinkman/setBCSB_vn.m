function S = setBCSB_vn(S,BC)
% setBCSB_vn - Implements the normal velocity boundary conditions into the system S
%
% SYNOPSIS:
%   S = setBCSB_vn(S,BC) 
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
  
  % Choose the normal velocity boundary conditions
  isVeln  = BC.nodetype == 1; %strcmp('velocity_n',BC.nodetype);
  BCdofs  = BC.node(isVeln);
  BCval   = BC.nodevalue(isVeln);
  BCtag   = BC.nodetag(isVeln);
  
  % Separate the upper, lower, left and right boundaries
  ile = find(BCtag==1); [BCdofs_le,i_le] = unique(BCdofs(ile));
  ir  = find(BCtag==2); [BCdofs_r,i_r]   = unique(BCdofs(ir));
  il  = find(BCtag==3); [BCdofs_l,i_l]   = unique(BCdofs(il));
  iu  = find(BCtag==4); [BCdofs_u,i_u]   = unique(BCdofs(iu));
  ib  = find(BCtag==5); [BCdofs_b,i_b]   = unique(BCdofs(ib));
  it  = find(BCtag==6); [BCdofs_t,i_t]   = unique(BCdofs(it));
  mb2 = mean(mean(abs(S.B2)));

  % The left boundary
  if(~isempty(ile))
    % move the respective columns to the right-hand side
    S.F1 = S.F1-S.B1(:,BCdofs_le)*BCval(i_le);
    if(isfield(S,'D1'))
      S.Q  = S.Q - S.D1(:,BCdofs_le)*BCval(i_le);
      S.D1(:,BCdofs_le) = 0;
    else
      S.Q  = S.Q-S.C1(BCdofs_le,:)'*BCval(i_le);
    end
    S.B1(:,BCdofs_le) = 0;
    S.B1(BCdofs_le,:) = 0;
    S.C1(BCdofs_le,:) = 0;
    % insert the velocity values into F 
    S.F1(BCdofs_le) = BCval(i_le)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_le)
      S.B1(BCdofs_le(i),BCdofs_le(i)) = mb2;
    end
  end;
  
  % The right boundary
  if(~isempty(ir)) 
    % move the respective columns to the right-hand side
    S.F1 = S.F1-S.B1(:,BCdofs_r)*BCval(i_r);
    if(isfield(S,'D1'))
      S.Q  = S.Q - S.D1(:,BCdofs_r)*BCval(i_r);
      S.D1(:,BCdofs_r) = 0;
    else    
      S.Q  = S.Q-S.C1(BCdofs_r,:)'*BCval(i_r);
    end
    S.B1(:,BCdofs_r) = 0;
    S.B1(BCdofs_r,:) = 0;
    S.C1(BCdofs_r,:) = 0;
    % insert the velocity values into F 
    S.F1(BCdofs_r) = BCval(i_r)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_r)
      S.B1(BCdofs_r(i),BCdofs_r(i)) = mb2;
    end
  end;

  % The lower boundary
  if(~isempty(il))
    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_l)*BCval(i_l);
    if(isfield(S,'D2'))
      S.Q  = S.Q-S.D2(:,BCdofs_l)*BCval(i_l);
      S.D2(:,BCdofs_l) = 0;
    else
      S.Q  = S.Q-S.C2(BCdofs_l,:)'*BCval(i_l);
    end
    S.B2(:,BCdofs_l) = 0;
    S.B2(BCdofs_l,:) = 0;
    S.C2(BCdofs_l,:) = 0;
    % insert the velocity values into F 
    S.F2(BCdofs_l) = BCval(i_l)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_l)
      S.B2(BCdofs_l(i),BCdofs_l(i)) = mb2;
    end
  end;
  
  % The upper boundary
  if(~isempty(iu))
    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_u)*BCval(i_u);
    if(isfield(S,'D2'))
      S.Q = S.Q-S.D2(:,BCdofs_u)*BCval(i_u);
      S.D2(:,BCdofs_u) = 0;
    else
      S.Q = S.Q-S.C2(BCdofs_u,:)'*BCval(i_u);
    end
    S.B2(:,BCdofs_u) = 0;
    S.B2(BCdofs_u,:) = 0;
    S.C2(BCdofs_u,:) = 0;
    % insert the velocity values into F 
    S.F2(BCdofs_u) = BCval(i_u)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_u)
      S.B2(BCdofs_u(i),BCdofs_u(i)) = mb2;
    end
  end;

  % The bottom boundary
  if(~isempty(ib))
    % move the respective columns to the right-hand side
    S.F3 = S.F3-S.B3(:,BCdofs_b)*BCval(i_b);
    if(isfield(S,'D3'))
      S.Q  = S.Q-S.D3(:,BCdofs_b)*BCval(i_b);
      S.D3(:,BCdofs_b) = 0;
    else
      S.Q  = S.Q-S.C3(BCdofs_b,:)'*BCval(i_b);
    end
    S.B3(:,BCdofs_b) = 0;
    S.B3(BCdofs_b,:) = 0;
    S.C3(BCdofs_b,:) = 0;
    % insert the velocity values into F 
    S.F3(BCdofs_b) = BCval(i_b)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_b)
      S.B3(BCdofs_b(i),BCdofs_b(i)) = mb2;
    end
  end

  % The top boundary
  if(~isempty(it))
    % move the respective columns to the right-hand side
    S.F3 = S.F3-S.B3(:,BCdofs_t)*BCval(i_t);
    if(isfield(S,'D3'))
      S.Q  = S.Q-S.D3(:,BCdofs_t)*BCval(i_t);
      S.D3(:,BCdofs_t) = 0;
    else  
      S.Q  = S.Q-S.C3(BCdofs_t,:)'*BCval(i_t);
    end
    S.B3(:,BCdofs_t) = 0;
    S.B3(BCdofs_t,:) = 0;
    S.C3(BCdofs_t,:) = 0;
    % insert the velocity values into F 
    S.F3(BCdofs_t) = BCval(i_t)*mb2;
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_t)
      S.B3(BCdofs_t(i),BCdofs_t(i)) = mb2;
    end
  end;
