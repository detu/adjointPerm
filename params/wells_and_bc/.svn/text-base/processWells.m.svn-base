function W = processWells(G, rock, grdecl, varargin)
% processWells -- Process output from function 'readGRDECL' to insert wells
% into the grid.
%
% SYNOPSIS:
%   W = processWells(G, rock,grdecl, varargin)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   rock    - Rock data structure.  Must contain valid field 'perm'.
%
%   grdecl  - Output from function 'readGRDECL', with field 'wells'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - InnerProduct -- The inner product with which to define
%                                 the mass matrix.
%                                 String.  Default value = 'ip_simple'.
%                 Supported values are:
%                   - 'ip_simple'
%                   - 'ip_tpf'
%                   - 'ip_quasitpf'
%                   - 'ip_rt'
% RETURNS:
%   W       - Updated (or freshly created) well structure, each element
%             of which has the following fields:
%               cells -- Grid cells perforated by this well.
%               type  -- Well control type.
%               val   -- Target control value.
%               r     -- Well bore radius.
%               dir   -- Well direction.
%               WI    -- Well index.
%               dZ    -- Vertical displacement of each well perforation
%                        measured from 'highest' horizontal contact (i.e.
%                        the 'TOP' contact with the minimum 'Z' value
%                        counted amongst all cells perforated by this
%                        well).
%               name  -- Well name.
%               compi -- Fluid composition--only used for injectors.
%
% SEE ALSO:
%   verticalWell, addWell, readGRDECL, readWellKW, processGRDECL.
%
% REMARKS:
%  * Sorting of well cells only works for structured grids/no gravity.
%
%  * Due to the way the well index (WI) is calculated in 'addWell' we
%    cannot allow varying radius and changing well directions when WI is
%    not supplied. (See get_diam_dir)
%
%  * When data for one completion in a well is given more than once, the
%    data first read is used. (See comp_cellInx)

%{
#COPYRIGHT#
%}

% $Id: processWells.m 2986 2009-10-13 10:25:47Z hnil $

   opt = struct('InnerProduct', 'ip_tpf');
   opt = merge_options(opt, varargin{:});

   ip = opt.InnerProduct;

   W = struct([]);

   % Loop through all wells
   for i = 1 : size(grdecl.wells,2),

      w = grdecl.wells(i);

      % Check if well is open, if not go to next well.
      if strcmp(w.state, 'SHUT')
         fprintf('Well %name is %s and will not be added to W', ...
                 name, w.state);
         continue
      end

      %% Define well parameters. ------------------------------------------
      %
      name = w.name;
      isInj = strcmp(w.type, 'Injection');

      % Composition of inj wells (ignored for prod wells).
      comp_i = [0, 0, 0];
      if isInj, comp_i = get_comp(w.inj_type); end

      % Well control and value.
      [type, val] = get_control(w, isInj);
      %convert units to SI
      if(strcmp(type,'bhp'))
          val = val*barsa;
      elseif(strcmp(type,'rate'))
          val = val/day;
      else
          error('Unknown control mode')
      end

      % Find active well cells and index to vector values of active cells.
      [cellInx, inx] = comp_cellInx(G, w.I, w.J, w.I_vec, w.J_vec, ...
                                    w.K_vec_u, w.K_vec_l, w.state_vec);

      %% Construct well. --------------------------------------------------
      %
      s = logical(w.state_vec);
      transm = w.transm_vec(s)*(0.00852701472/(day*barsa));

      % CASE 1: transmissibility (WI) is defined for all cells:
      if all(transm > 0),
        if(~isempty(cellInx))
          W = addWell(G, rock, W, cellInx, 'Type', type, 'Val', val,     ...
                      'Name', name, 'Comp_i', comp_i, 'WI', transm(inx), ...
                      'InnerProduct', ip,'refDepth',min(G.cells.centroids(cellInx,3)));


        else
          display(['Warning: skipping well without perforation in domain']); 
        end
         
      else
      % CASE 2: transmissibility (WI) is not defined for all cells.. :
      % For now : only allow same diameter in entire well and and either
      % horizontal or vertical well.
      % units???
         Kh = w.Kh_vec(s)*milli*darcy;
         skin = w.skin_vec(s);

         [diam, dir] = get_diam_dir(w);
         if(~isempty(cellInx))
            W = addWell(G, rock, W, cellInx, ...
                         'InnerProduct', ip, 'Type', type, 'Val', val, ...
                         'Radius', diam(inx)/2, 'Dir', dir(inx), 'Name', name,   ...
                         'Comp_i', comp_i, 'WI', transm(inx), 'Kh',    ...
                        Kh(inx), 'Skin', skin(inx));
         else
          display(['Warning: skipping well without perforation in domain']); 
        end
      end
   end
end


%--------------------------------------------------------------------------
% Private helper functions follow
%--------------------------------------------------------------------------


function comp = get_comp(inj_type)
   switch lower(inj_type)
      case 'water'
         comp = [1, 0 , 0];
      case 'oil'
         comp = [0, 1 , 0];
      case 'gas'
         comp = [0, 0 , 1];
   end
end

%--------------------------------------------------------------------------

function [type, val] = get_control(w, isInj)
   switch upper(w.control)
      case {'BHP'}
         type = 'bhp';
         val = w.bhp;
      case {'RESV'} % reservoir fluid vol_rate
         type = 'rate';
         if isInj,
            val  = w.vol_rate;
         else
            val  = -w.vol_rate;
         end
         
      case {'RATE'}, % reservoir fluid vol_rate
       error('RATE is not supported in the solvers');
         type = 'rate';
         assert(isInj)
         val  = w.surf_vol_rate;  
               
     case {'LRAT'}, % reservoir fluid vol_rate
         type = 'rate';
         assert(~isInj)
         val  = -w.liquid_rate;

           
      otherwise
        error('processWells:ControlMode:NotSupported', ...
              'Control mode ''%s'' is not supported\n', w.control);
   end
end

%--------------------------------------------------------------------------

function [diam, dir] = get_diam_dir(w)
   % For now: check well diameter (only allow same diameter in entire well)
   % assert that the radius is the same for all completions (all COMPDAT).
   num_d = numel(w.diam_vec);
   %diam = w.diam_vec(1);
   diam = w.diam_vec';
   dir = w.dir_vec';
   %assert(diam*num_d == sum(w.diam_vec))
   %assert(abs(diam*num_d -sum(w.diam_vec))/sum(w.diam_vec) < 1e-4)

   % For now: only allow either horizontal or vertical well
   % i.e., 'dir_vec' only holds one type of letter.
%    num = sum(1 : numel(w.dir_vec));
% 
%    if     sum(regexp(w.dir_vec, 'Z')) == num,
%       dir = 'Z';
%    elseif sum(regexp(w.dir_vec, 'Y')) == num,
%       dir = 'Y';
%    elseif sum(regexp(w.dir_vec, 'X')) == num,
%       dir = 'X';
%    else
%       % not really needed.. is enough that the direction is the same for
%       % those that must be calculated?
%      error('processWells:wellDirection:NotSupported',                      ...
%           ['The well direction must be the same for all COMPDAT\n',        ...
%            'when transmissibilites are not given for all completions\n. ', ...
%            'Combination %s is not supported\n'], w.dir_vec);
%    end
end
%--------------------------------------------------------------------------

function [cells, inx] = comp_cellInx(G, I, J, I_vec, J_vec, K_u, K_l, state)
   % NB: must consider state when using inx!
   % NB: return index (inx) to valid compdats in case one completion is
   % specified in several compdats..
   assert (isfield(G      , 'cartDims') && (numel(G.cartDims) == 3));
   assert (isfield(G.cells, 'indexMap'));

   inx   = []; %inx to use in transm_vec, Kh_vec or skin_vec
   cells = [];

   s = logical(state);
   [I_vec, J_vec, K_u, K_l] = deal(I_vec(s), J_vec(s), K_u(s), K_l(s));
   I_vec(I_vec <0) = I;
   J_vec(J_vec <0) = J;


   if (min(K_u) < 1) || (max(K_l) > G.cartDims(3)) || any(K_u>K_l),
      error(msgid('CompSpec:Invalid'), ...
            'Vertical completion specified outside model.');
   end
   if any(I_vec < 1) || any(I_vec > G.cartDims(1)) || ...
      any(J_vec < 1) || any(J_vec > G.cartDims(2)),
      error(msgid('CompSpec:Invalid'), ...
            'Horizontal completion specified outside model.');
   end

   for i = 1 : numel(I_vec),
      I = I_vec(i);
      J = J_vec(i);
      K = K_u(i) : K_l(i);
      num = numel(1);

      % Taken from verticalWell() -----
      assert (all([numel(I), numel(J)] == 1));

      nl = numel(K); I = I(ones([nl, 1])); J = J(ones([nl, 1]));
      ix = sub2ind(reshape(G.cartDims,1,[]), I, J, K(:));

      wc     = false([prod(G.cartDims), 1]);
      wc(ix) = true;
      wc     = find(wc(G.cells.indexMap));
      % -----
      cells = [cells; wc];
      if(~isempty(wc))
        inx = [inx; i*ones(num,1)];
      end
   end

   % Cells are sorted.
   [cells, first_inx, last_inx] = unique(cells, 'first');

   if first_inx ~= last_inx
      warning('processWells:comp_cellInx', ...
              'Ignoring several compdats for same completion.')
      inx = inx(first_inx);
   end
end
