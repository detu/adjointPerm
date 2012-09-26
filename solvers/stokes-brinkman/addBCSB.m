function bc = addBCSB(bc, f, t, v, G, Dofs, varargin)
% addBCSB - Add (Stokes-Brinkman) boundary conditions to (new or existing) BC object
%
% SYNOPSIS:
%   bc = addBCSB(bc, faces, type, values, G, Dofs)
%   bc = addBCSB(bc, faces, type, values, G, Dofs, 'pn', pv, ...)
%
% PARAMETERS:
%   bc      - Boundary condition structure from a prior call to 'addBCSB'
%             which will be updated on output or an empty array (bc==[]) in
%             which case a new boundary condition structure is created.
%
%   faces   - Global faces in external model for which this boundary
%             condition should be applied.
%
%   type    - Type of boundary condition.  Supported values are 'pressure',
%             'velocity_n' and 'velocity_t'.
%
%   values  - Boundary condition value.  Interpreted as a pressure value
%             when type=='pressure', as a normal velocity value when 
%             type=='velocity_n' and as a tangential velocity value when
%             type=='velocity_t'. One scalar value for each node in 'nodes'.
%
%   G       - Grid structure
%
%   Dofs    - Degrees-of-freedom structure  
%
%   'pn'/pv - List of propertyname/propertyvalue pairs.
%             PropertyNames    PropertyValues     
%             'facetag'        The face tag faces in 'faces'
%             'node'           Nodes connected to 'faces'
%             'nodevalue'      BC value of nodes connected to 'faces' 
%             'nodetag'        The tag for nodes for which BCs are procided     
%
% RETURNS:
%   bc - New or updated boundary condition structure having the following
%        fields:
%           - face       -- External faces for which explicit BCs are provided.
%           - type       -- Array of int denoting type of BC: 
%                           1 = velocity_n, 2 = pressure, 3 = velocity_t.
%           - value      -- Boundary condition values for all faces in 'face'.
%           - node       -- External nodes for which explicit BCs are provided.
%           - nodevalue  -- Boundary condition values for all nodes in 'node'.
%           - nodetag    -- Tag for external nodes for which explicit BCs
%                           are provided. 
%           - nodetype   -- Array of int denoting type of BC for nodes 
%                           (see 'type' above).
%
  
error(nargchk(6, 14, nargin, 'struct'));

if isempty(bc),
    bc = struct('face', [], 'type', [], 'value', [],...
                'node', [], 'nodetype', [], 'nodevalue', [], 'nodetag', []);                
end
opt = struct('facetag',   [], ...
             'node',      [], ...
             'nodevalue', [], ...
             'nodetag',   [], ...
             'nodetype',  []);
opt = merge_options(opt, varargin{:});

assert (numel(f) == numel(v));

dim    = numel(G.cartDims);
dim_no = size(Dofs.face2nodes,2);


if strcmp(t, 'velocity_n')
   t = 1;
elseif strcmp(t,  'pressure')
   t = 2;
else % velocity_t
   t = 3;
end
[f,ind] = sort(f);
v       = v(ind);

if isempty(opt.facetag)
  dum     = false(G.faces.num,1);   
  dum(f)  = true;    
  tags    = G.cellFaces(dum(G.cellFaces(:,1)),2);
  f       = G.cellFaces(dum(G.cellFaces(:,1)),1);
  [f,ind] = sort(f);
  tags    = tags(ind);

  if isempty(opt.nodetag)
    ntag = reshape(repmat(tags,1,dim_no)',[],1);
  else
    ntag = opt.nodetag;
  end

  if isempty(opt.node)
    n = reshape(Dofs.face2nodes(f,:)',[],1);
  else
    n = opt.node;
  end

  if isempty(opt.nodevalue)
    nv = reshape(repmat(v,1,dim_no)',[],1);
  else
    nv = opt.nodevalue;
  end
else
  opt.facetag = opt.facetag(ind);
  
  if isempty(opt.nodetag)
    ntag = reshape(repmat(opt.facetag,1,dim_no)',[],1);
  else
    ntag = opt.nodetag;
  end
    
  if isempty(opt.node)
    n = reshape(Dofs.face2nodes(f,:)',[],1);
  else
    n = opt.node;
  end
    
  if isempty(opt.nodevalue)
    nv = reshape(repmat(v,1,dim_no)',[],1);
  else
    nv = opt.nodevalue;
  end

end

bc.face      = [bc.face       ; f];
bc.type      = [bc.type       , t*ones([1, numel(f)])];
bc.value     = [bc.value      ; v];
bc.node      = [bc.node       ; n];
bc.nodevalue = [bc.nodevalue  ; nv];
bc.nodetag   = [bc.nodetag    ; ntag];
bc.nodetype  = [bc.nodetype, t*ones([1, numel(n)])];
