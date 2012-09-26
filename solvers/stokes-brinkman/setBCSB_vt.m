function S = setBCSB_vt(S,BC)
% setBCSB_vt - Implements the tangential velocity boundary conditions into the system S
%
% SYNOPSIS:
%   S = setBCSB_vt(S,BC) 
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

  % Choose the tangential velocity boundary conditions
  isVelt = strcmp('velocity_t',BC.nodetype);
  BCdofs = BC.node(isVelt);
  BCval  = BC.nodevalue(isVelt);
  BCtag  = BC.nodetag(isVelt);

  if isfield(S, 'B3')
    dim=3;
  else
    dim=2;
  end
    
  % Separate the upper, lower, left and right boundaries
  ile = find(BCtag==1); [BCdofs_le,i_le] = unique(BCdofs(ile));
  ir  = find(BCtag==2); [BCdofs_r,i_r]   = unique(BCdofs(ir));
  il  = find(BCtag==3); [BCdofs_l,i_l]   = unique(BCdofs(il));
  iu  = find(BCtag==4); [BCdofs_u,i_u]   = unique(BCdofs(iu));
  ib  = find(BCtag==5); [BCdofs_b,i_b]   = unique(BCdofs(ib));
  it  = find(BCtag==6); [BCdofs_t,i_t]   = unique(BCdofs(it));
  mb2 = mean(mean(abs(S.B2))); 

  % The left boundary - y and z tangential
  if(~isempty(BCdofs_le))

    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_le)*BCval(i_le);
    if(isfield(S,'D2'))
      S.Q  = S.Q-S.D2(:,BCdofs_le)*BCval(i_le);
      S.D2(:,BCdofs_le) = 0;
    else   
      S.Q  = S.Q-S.C2(BCdofs_le,:)'*BCval(i_le);
    end
    S.B2(:,BCdofs_le) = 0;
    S.B2(BCdofs_le,:) = 0;
    S.C2(BCdofs_le,:) = 0;
    % insert the velocity values into F 
    S.F2(BCdofs_le) = BCval(i_le)*mb2;
    
    if dim==3
      % move the respective columns to the right-hand side
      S.F3 = S.F3-S.B3(:,BCdofs_le)*BCval(i_le);
      if(isfield(S,'D3'))
        S.Q  = S.Q-S.D3(:,BCdofs_le)*BCval(i_le);
        S.D3(:,BCdofs_le) = 0;
      else
        S.Q  = S.Q-S.C3(BCdofs_le,:)'*BCval(i_le);
      end
      S.B3(:,BCdofs_le) = 0;
      S.B3(BCdofs_le,:) = 0;
      S.C3(BCdofs_le,:) = 0;
      % insert the velocity values into F 
      S.F3(BCdofs_le) = BCval(i_le)*mb2;
    end
    
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_le)
      S.B2(BCdofs_le(i),BCdofs_le(i)) = mb2;
      if dim==3
        S.B3(BCdofs_le(i),BCdofs_le(i)) = mb2;
      end
    end
    
  end;
  
  % The right boundary - y and z tangential
  if(~isempty(BCdofs_r))
    
    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_r)*BCval(i_r);
    if(isfield(S,'D2'))
      S.Q  = S.Q-S.D2(:,BCdofs_r)*BCval(i_r);
      S.D2(:,BCdofs_r) = 0;
    else
      S.Q  = S.Q-S.C2(BCdofs_r,:)'*BCval(i_r);
    end
    S.B2(:,BCdofs_r) = 0;
    S.B2(BCdofs_r,:) = 0;
    S.C2(BCdofs_r,:) = 0;
    % insert the velocity values into F 
    S.F2(BCdofs_r) = BCval(i_r)*mb2;
    
    if dim==3
      % move the respective columns to the right-hand side
      S.F3 = S.F3-S.B3(:,BCdofs_r)*BCval(i_r);
    if(isfield(S,'D3'))
      S.Q  = S.Q-S.D3(:,BCdofs_r)*BCval(i_r);
      S.D3(:,BCdofs_r) = 0;
    else  
      S.Q  = S.Q-S.C3(BCdofs_r,:)'*BCval(i_r);
    end
      S.B3(:,BCdofs_r) = 0;
      S.B3(BCdofs_r,:) = 0;
      S.C3(BCdofs_r,:) = 0;
      % insert the velocity values into F 
      S.F3(BCdofs_r) = BCval(i_r)*mb2;
    end  
    
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_r)
      S.B2(BCdofs_r(i),BCdofs_r(i)) = mb2;
      if dim==3
        S.B3(BCdofs_r(i),BCdofs_r(i)) = mb2;
      end
    end
    
  end;
  
  % The lower boundary - x and z tangential
  if(~isempty(BCdofs_l))
    % move the respective columns to the right-hand side
    S.F1 = S.F1-S.B1(:,BCdofs_l)*BCval(i_l);
    if(isfield(S,'D1'))
      S.Q  = S.Q-S.D1(:,BCdofs_l)*BCval(i_l);
      S.D1(:,BCdofs_l) = 0;
    else   
      S.Q  = S.Q-S.C1(BCdofs_l,:)'*BCval(i_l);
    end
    S.B1(:,BCdofs_l) = 0;
    S.B1(BCdofs_l,:) = 0;
    S.C1(BCdofs_l,:) = 0;
    % insert the velocity values into F 
    S.F1(BCdofs_l) = Cval(i_l)*mb2;
    
    if dim==3
      % move the respective columns to the right-hand side
      S.F3 = S.F3-S.B3(:,BCdofs_l)*BCval(i_l);
      if(isfield(S,'D3'))
        S.Q  = S.Q-S.D3(:,BCdofs_l)*BCval(i_l);
        S.D3(:,BCdofs_l) = 0;
      else 
        S.Q  = S.Q-S.C3(BCdofs_l,:)'*BCval(i_l);
      end
      S.B3(:,BCdofs_l) = 0;
      S.B3(BCdofs_l,:) = 0;
      S.C3(BCdofs_l,:) = 0;
      % insert the velocity values into F 
      S.F3(BCdofs_l) = BCval(i_l)*mb2;
    end
    
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_l)
      S.B1(BCdofs_l(i),BCdofs_l(i)) = mb2;
      if dim==3
        S.B3(BCdofs_l(i),BCdofs_l(i)) = mb2;
      end
    end
    
  end;

  % The upper boundary - x and z tangential
  if(~isempty(BCdofs_u))
   % move the respective columns to the right-hand side
   S.F1 = S.F1-S.B1(:,BCdofs_u)*BCval(i_u);
   if(isfield(S,'D1'))
     S.Q  = S.Q-S.D1(:,BCdofs_u)*BCval(i_u);
     S.D1(:,BCdofs_u) = 0;
   else
     S.Q  = S.Q-S.C1(BCdofs_u,:)'*BCval(i_u);
   end
   S.B1(:,BCdofs_u) = 0;
   S.B1(BCdofs_u,:) = 0;
   S.C1(BCdofs_u,:) = 0;
   % insert the velocity values into F 
   S.F1(BCdofs_u) = BCval(i_u)*mb2;
    
   if dim==3
     % move the respective columns to the right-hand side
     S.F3 = S.F3-S.B3(:,BCdofs_u)*BCval(i_u);
     if(isfield(S,'D3'))
       S.Q  = S.Q-S.D3(:,BCdofs_u)*BCval(i_u);
       S.D3(:,BCdofs_u) = 0;
     else   
       S.Q  = S.Q-S.C3(BCdofs_u,:)'*BCval(i_u);
     end
     S.B3(:,BCdofs_u) = 0;
     S.B3(BCdofs_u,:) = 0;
     S.C3(BCdofs_u,:) = 0;
     % insert the velocity values into F 
     S.F3(BCdofs_u) = BCval(i_u)*mb2;
   end
   
   % set B(u_k,u_k) equal to 1
   for i = 1:numel(BCdofs_u)
     S.B1(BCdofs_u(i),BCdofs_u(i)) = mb2;
     if dim==3
       S.B3(BCdofs_u(i),BCdofs_u(i)) = mb2;
     end
   end
   
  end;
  
  % The bottom boundary - x and y tangential
  if(~isempty(BCdofs_b))
    
    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_b)*BCval(i_b);
    S.F1 = S.F1-S.B1(:,BCdofs_b)*BCval(i_b);
    if(isfield(S,'D1'&isfield(S,'D2')))
      S.Q  = S.Q-S.D2(:,BCdofs_b)*BCval(i_b);
      S.Q  = S.Q-S.D1(:,BCdofs_b)*BCval(i_b);
      S.D2(BCdofs_b,:) = 0;
      S.D1(BCdofs_b,:) = 0;
    else  
      S.Q  = S.Q-S.C2(BCdofs_b,:)'*BCval(i_b);
      S.Q  = S.Q-S.C1(BCdofs_b,:)'*BCval(i_b);
    end
    S.B2(:,BCdofs_b) = 0;
    S.B1(:,BCdofs_b) = 0;
    S.B2(BCdofs_b,:) = 0;
    S.B1(BCdofs_b,:) = 0;
    S.C2(BCdofs_b,:) = 0;
    S.C1(BCdofs_b,:) = 0;
    % insert the velocity values into F 
    S.F2(BCdofs_b) = BCval(i_b)*mb2;
    S.F1(BCdofs_b) = BCval(i_b)*mb2;
   
    % set B(u_k,u_k) equal to 1
    for i = 1:numel(BCdofs_b)
      S.B2(BCdofs_b(i),BCdofs_b(i)) = mb2;
      S.B1(BCdofs_b(i),BCdofs_b(i)) = mb2;
    end
    
  end;
  
  % The top boundary - x and y tangential
  if(~isempty(BCdofs_t))

    % move the respective columns to the right-hand side
    S.F2 = S.F2-S.B2(:,BCdofs_t)*BCval(i_t);
    S.F1 = S.F1-S.B1(:,BCdofs_t)*BCval(i_t);
     if(isfield(S,'D1'&isfield(S,'D2')))
       S.Q  = S.Q-S.D2(:,BCdofs_t)*BCval(i_t);
       S.Q  = S.Q-S.D1(:,BCdofs_t)*BCval(i_t);
       S.D2(BCdofs_t,:) = 0;
       S.D1(BCdofs_t,:) = 0;
     else
       S.Q  = S.Q-S.C2(BCdofs_t,:)'*BCval(i_t);
       S.Q  = S.Q-S.C1(BCdofs_t,:)'*BCval(i_t);
     end
     S.B2(:,BCdofs_t) = 0;
     S.B1(:,BCdofs_t) = 0;
     S.B2(BCdofs_t,:) = 0;
     S.B1(BCdofs_t,:) = 0;
     S.C2(BCdofs_t,:) = 0;
     S.C1(BCdofs_t,:) = 0;
     
     % insert the velocity values into F 
     S.F2(BCdofs_t) = BCval(i_t)*mb2;
     S.F1(BCdofs_t) = BCval(i_t)*mb2;
     
     % set B(u_k,u_k) equal to 1
     for i = 1:numel(BCdofs_t)
       S.B2(BCdofs_t(i),BCdofs_t(i)) = mb2;
       S.B1(BCdofs_t(i),BCdofs_t(i)) = mb2;
     end
     
  end;
  