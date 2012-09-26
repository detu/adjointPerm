function p2 = refinePartitionForWells(G, W, p, varargin)
% refinePartitionForWells -- Refine coarse block partitions near wells
%
% SYNOPSIS:
%   partition = refinePartitionForWells(G, W, partition)
%   partition = refinePartitionForWells(G, W, partition, 'pn', pv, ...)
%
% PARAMETERS:
%   G         - grid_structure structure.  Currently assumed to stem from a
%               logically Cartesian grid (i.e. must either have a field
%               'G.cells.ijkMap' or a field 'G.cartDims').
%
%   W         - Well structure as built by addWell &c.
%
%   partition - Existing partition, created e.g. by function `partitionUI'.
%
%   'pn'/pv   - List of 'key'/value pairs for supplying optional
%               parameters.  All optional parameters go towards specifying
%               a vector `r' of refinement layer radii of the form
%
%                   r = [r0, r1, r2, ..., rf]     % (*)
%
%               All radii, except possibly r0, must be positive.  The
%               radius r0 of the inner-most refinement layer about a well
%               must be non-negative.  If r0==0 each well in W is placed
%               into separate coarse blocks.
%
%               The supported options and values are:
%                 - explRef -- Explicit setting of layer radii, i.e. user
%                              provides vector 'r' of (*) directly.
%                              Overrides all other options.
%
%                 - nLayer  -- Number of "automatic" refinement layers
%                              around each well.  Default value: nLayer=4.
%                              Setting nLayer=0 prevents near-well
%                              refinement. Setting nLayer=1 means that each
%                              well is placed into its own coarse block,
%                              but no other refinement is attempted.
%
%                 - refStrat-- Refinement strategy.  String.  Which named,
%                              predefined refinement strategy should be
%                              employed.  Default: refStrat = 'geometric'.
%
%                              The currently known strategies are:
%                                - 'geometric':
%                                     Geometric (power-of-two) increases in
%                                     refinement layer radius.  The default
%                                     of nLayer=4 yields r = [0, 1, 2, 4].
%                                - 'fibonacci':
%                                     Increases refinement layer radius
%                                     according to the well-known Fibonacci
%                                     sequence.  The default of nLayer=4
%                                     yields r = [0, 1, 1, 2].
%
% EXAMPLE:
%   G = buildMatchingGrid(grdecl);
%   p = partitionUI(G, coarseDim);
%   p = refinePartitionForWells(G, W, p, 'explRef', [0, pow2(1 : 3)]);
%   p = processPartition(G, p);
%
% RETURNS:
%   partition - Refined partition.
%
% SEE ALSO:
%   partitionCartGrid, partitionUI, processPartition, generateCoarseGrid.

% $Id: refinePartitionForWells.m 1773 2009-03-19 12:47:11Z bska $

%% Check arguments
isCart = isstruct(G) && ...
         ((isfield(G, 'cells') && ...
           (isfield(G.cells, 'ijkMap') || ...
            isfield(G.cells, 'indexMap'))) || ...
          isfield(G, 'cartDims'));
if ~isCart,
   error(id('GRID:NOIJK'), 'Grid must be logically Cartesian');
end

%% Setup
p2 = p;
wc = { W.cells };
nw = numel(W);
nb = max(p);

ijkMap = getCellIJK(G);

opt = struct('nLayer', 4, 'refStrat', 'geometric', 'explRef', []);
opt = merge_options(opt, varargin{:});
r   = defineRefinement(opt);

%% Refinement loop
if ~isempty(r),
   for w = 1 : nw,
      switch lower(W(w).dir),
         case 'z', [d0, d1, d2] = deal(3, 1, 2);
         case 'x', [d0, d1, d2] = deal(1, 2, 3);
         case 'y', [d0, d1, d2] = deal(2, 3, 1);
         otherwise,
            error(id('Well:Direction'),                       ...
                  'Unknown well direction ''%s'' in well %d', ...
                  W(w).dir, w);
      end

      wb   = p(wc{w});
      uwb  = unique(wb);
      nwb  = numel(uwb);
      wijk = double(ijkMap(wc{w}, :));
      if norm(wijk(:,[d1,d2]) - ...
              repmat(wijk(1, [d1,d2]), [numel(wb), 1]), ...
              inf) > 0,
         error(id('Well:IJK'), ...
               'Unexpected non-constant IJK for well %d', w);
      end

      for b = 1 : nwb,
         inblock = p(wc{w}) == uwb(b);
         layers = wijk(inblock, d0);
         k     = 0; old_cells = [];
         cells = wc{w}(inblock);
         for j = 1 : numel(r),
            while k < r(j),
               cells = neighboursByNodes(G, cells);
               k = k + 1;
            end
            % Exclude cells outside the required 'd0'-layers
            cells = cells(ismember(ijkMap(cells,d0), layers));

            nb = nb + 1;
            p2(setdiff(cells, old_cells)) = nb;

            old_cells = cells;
         end
      end
   end

   %% Remove empty coarse blocks if any
   p2 = compressPartition(p2);
end


%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------


function r = defineRefinement(opt)
r = [];
if ~isempty(opt.explRef),
   r = opt.explRef;
elseif opt.nLayer > 0,
   n = opt.nLayer;
   switch lower(opt.refStrat),
      case {'geometric', 'pow2'},
         r = [0, pow2((1 : n - 1) - 1)];
      case {'fibonacci'},
         r = zeros([1, n]);
         if n > 1,
            r(2) = 1;
            for k = 3 : n, r(k) = r(k-2) + r(k-1); end
         end
      otherwise
         error(id('Strategy:Unsupported'),              ...
               ['Refinement strategy ''', opt.refStrat, ...
                ''' is not supported']);
   end
end

if ~isempty(r),
   if r(1) < 0,
      error(id('R0:Negative'), ...
            'Negative inner refinement radius is unsupported');
   elseif numel(r) > 1 && any(r(2:end) < 1),
      error(id('Ri:NonPositive'), ...
            'All non-inner refinement radii must be positive');
   end
end

r = cumsum(r);

%--------------------------------------------------------------------------

function s = id(s)
s = ['refinePartitionForWell:', s];

%--------------------------------------------------------------------------

function ijkMap = getCellIJK(G)
if isfield(G.cells, 'ijkMap'),
   ijkMap = G.cells.ijkMap;
else
   if isfield(G.cells, 'indexMap'),
      im = G.cells.indexMap;
   else
      im = (1 : G.cells.num) .';
   end
   [ijkMap{1:3}] = ind2sub(G.cartDims(:).', im);
   ijkMap = [ijkMap{:}];
end
