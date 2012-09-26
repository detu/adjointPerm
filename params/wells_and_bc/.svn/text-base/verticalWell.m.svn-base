function W = verticalWell(W, G, rock, I, J, K, varargin)
%Insert a vertical well into the simulation model.
%
% SYNOPSIS:
%   W = verticalWell(W, G, rock, I, J, K)
%   W = verticalWell(W, G, rock, I, J, K, 'pn', pv)
%
% PARAMETERS:
%   W       - Input well structure defining existing wells.  Pass an empty
%             structure if there are no previously defined wells.  This
%             structure is updated on output.
%
%   G       - Grid data structure.  Must contain valid field
%             'G.cells.centroids'.  Call function 'computeGeometry' to
%             obtain these values.
%
%   rock    - Rock data structure.  Must contain valid field 'rock.perm'.
%
%   I, J    - Horizontal location of well head, measured in logical cell
%             indices.  'I' is the index along the first logical direction
%             while 'J' is the index along the second logical direction.
%
%   K       - A vector of layers in which this well should be completed.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - InnerProduct -- The inner product with which to define
%                                 the mass matrix. String.
%                                 Default value = 'ip_simple'.
%                 Supported values are:
%                   - 'ip_simple'
%                   - 'ip_tpf'
%                   - 'ip_quasitpf'
%                   - 'ip_rt'
%
%               - Type        -- Well control type.
%                                String.  Default value is 'bhp'.
%                 Supported values are:
%                   - 'bhp'  - Well is controlled by bottom hole pressure
%                              target.
%                   - 'rate' - Well is controlled by total rate target.
%
%               - Val         -- Pressure or rate value. Default value is 0
%
%               - Radius      -- Well bore radius [in meters].
%                                Default value is r = 0.1
%
%               - Name        -- Well name (string). OPTIONAL.
%
%               - Comp_i      -- Fluid composition for injection wells.
%                                Vector of saturations.
%                                Default value:
%                                  Comp_i = [1, 0, 0] (water injection)
%
%               - WI          -- Well productivity index.  Vector of length
%                                numC=numel(cellInx).
%                                Default value: WI = REPMAT(-1, [numC, 1]),
%                                whence the productivity index will be
%                                computed from available grid block data in
%                                grid blocks containing well completions.
%
%               - Kh          -- Permeability thickness. Vector of
%                                length numC=numel(cellInx).
%                                Default value: Kh = REPMAT(-1, [numC, 1]),
%                                whence the thickness will be computed from
%                                available grid block data in grid blocks
%                                containing well completions.
%
%               - Skin        -- Skin factor for computing effective well
%                                bore radius. Scalar value or vector of
%                                length numel(cellInx).
%                                Default value: 0.0 (no skin effect).
%
% RETURNS:
%   W - Updated well structure.
%
% EXAMPLE:
%   G = cartGrid([60, 220, 85], 0.3048.*[20, 10, 2].*[60, 220, 85]);
%   G = computeGeometry(G);
%   rock.perm  = repmat([500, 500, 50].*milli*darcy(), [G.cells.num, 1]);
%   W = struct([]);
%   W = verticalWell(W, G, rock,  1,   1, (1:85), ...
%                    'Type', 'bhp', 'Val', 300,   ...
%                    'Radius', 0.125, 'Name', 'P1');
%   W = verticalWell(W, G, rock, 60,   1, (1:85), ...
%                    'Type', 'bhp', 'Val', 300,   ...
%                    'Radius', 0.125, 'Name', 'P2');
%   W = verticalWell(W, G, rock,  1, 220, (1:85), ...
%                    'Type', 'bhp', 'Val', 300,   ...
%                    'Radius', 0.125, 'Name', 'P3');
%   W = verticalWell(W, G, rock, 60, 220, (1:85), ...
%                    'Type', 'bhp', 'Val', 300,   ...
%                    'Radius', 0.125, 'Name', 'P4');
%   W = verticalWell(W, G, rock, 30, 110, (1:85), ...
%                    'Type', 'bhp', 'Val', 300,   ...
%                    'Radius', 0.125, 'Name', 'I1', ...
%                    'Comp_i', [0.5, 0.5, 0]);
%
% SEE ALSO:
%   addWell, assembleWellSystem, addSource.

%{
#COPYRIGHT#
%}

% $Id: verticalWell.m 2338 2009-06-05 17:19:30Z bska $

assert (isfield(G      , 'cartDims') && (numel(G.cartDims) == 3));
assert (isfield(G.cells, 'indexMap'));
assert (all([numel(I), numel(J)] == 1));

if (min(K) < 1) || (max(K) > G.cartDims(3)),
   error(msgid('CompSpec:Invalid'), ...
         'Vertical completion specified outside model.');
end
if (I < 1) || (I > G.cartDims(1)) || (J < 1) || (J > G.cartDims(2)),
   error(msgid('CompSpec:Invalid'), ...
         'Horizontal completion specified outside model.');
end

error(nargchk(6, inf, nargin, 'struct'));


nl = numel(K); I = I(ones([nl, 1])); J = J(ones([nl, 1]));
ix = sub2ind(reshape(G.cartDims,1,[]), I, J, K(:));

wc     = false([prod(G.cartDims), 1]);
wc(ix) = true;
wc     = find(wc(G.cells.indexMap));

numC = numel(wc);
opt = struct('InnerProduct', 'ip_simple',                  ...
             'Dir'         , 'z',                          ...
             'Name'        , sprintf('W%d', numel(W) + 1), ...
             'Radius'      , 0.1,                         ...
             'Type'        , 'bhp',                        ...
             'Val'         , 0,                            ...
             'Comp_i'      , [1, 0, 0],                    ...
             'WI'          , repmat(-1, [numC, 1]),        ...
             'Kh'          , repmat(-1, [numC, 1]),        ...
             'Skin'        , zeros([numC, 1]));
opt = merge_options(opt, varargin{:});

W  = addWell(G, rock, W, wc, 'Type', opt.Type, 'Val', opt.Val,  ...
            'Radius', opt.Radius, 'Dir', 'z', 'Name', opt.Name, ...
            'InnerProduct', opt.InnerProduct, 'Comp_i', opt.Comp_i, ...
            'WI', opt.WI, 'Kh', opt.Kh, 'Skin', opt.Skin);
