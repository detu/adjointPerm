function relperm = initBlackOilRelPerm(T, varargin)
%Construct Black Oil relperm evaluation function.
%
% SYNOPSIS:
%   relperm = initBlackoilPVT(table)
%
% PARAMETERS:
%   table   - Black Oil PVT/Relperm table structure as defined by function
%             'readpvt'.
%
% RETURNS:
%   relperm - Function for evaluating relative permeability curves
%             specified by SWOF and SGOF.  Specifically, the call
%
%                [kr, dkr] = relperm(s)
%
%             evaluates the relative permeabilities (kr) and the
%             differentiated relative permeabilities (dkr) of the current
%             fluid composition (s, saturation).  The columns of 's' are
%             assumed to be
%
%                [Aqua, Liquid, Vapour]
%
%             and the return values follow the same convention.
%
% SEE ALSO:
%   initBlackoilFluid, initBlackOilPVT, swof, sgof.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
% $Revision: 2926 $

   opt = struct('verbose', false, 'phases', {{'A', 'L', 'V'}});
   opt = merge_options(opt, varargin{:});

   if ischar(opt.phases), opt.phases = num2cell(opt.phases); end

   A = find(strcmpi(opt.phases, 'A'));
   L = find(strcmpi(opt.phases, 'L'));
   V = find(strcmpi(opt.phases, 'V'));

   if numel(A)>1 || numel(L)>1 || numel(V)>1,
      error ('Cannot specify any phase more than once.');
   end
   illegal = ~ismember([opt.phases{:}], 'ALV');
   if any(illegal),
      s = 's'; if sum(illegal) == 1, s = ''; end
      error('Unknown phase%s specified', s);
   end

   if numel([A,L,V]) == 0, error('Huh!?, No phases?'); end

   % Get number of pvt regions
   num_pvt_regions = 1;
   if isfield(T, 'pvtnum'), num_pvt_regions = max(T.pvtnum);   end
   if num_pvt_regions > 1, error ('PVTNUM not supported yet'); end

   if isfield(T, 'satnum') && max(T.satnum(:)) > 1,
      error('Only a single saturation region is currently supported.');
   end

   Swco = 0;
   if isfield(T, 'swof'),
      [krw, krow, Swco, krocw] = swof(T);
   end

   if isfield(T, 'sgof'),
      [krg, krog] = sgof(T, Swco);
   end

   kro = @(s) default_kro(s, krow, krog, Swco, krocw);

   function [kr, varargout] = blackOilRelPerm(s)
      kr = zeros([size(s,1), 3]);

      [krw_s{1 : nargout}] = krw(s(:,  A       ));
      [kro_s{1 : nargout}] = kro(s(:, [A, L, V]));
      [krg_s{1 : nargout}] = krg(s(:,        V ));

      kr(:, [A, L, V]) = [krw_s{1}, kro_s{1}, krg_s{1}];

      if nargout > 1,
         dkr = zeros([size(s,1), 9]);
         col = @(i) (i - 1) * 3;

         dkr(:, col(A) +  A       ) = krw_s{2};
         dkr(:, col(L) + [A, L, V]) = kro_s{2};
         dkr(:, col(V) +        V ) = krg_s{2};

         varargout{1} = spblockdiag(dkr.', [3, 3]) .';
      end
   end

   relperm = @blackOilRelPerm;
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function S = spblockdiag(v, blocksize)
   % Construct sparse block diagonal matrix with (uniform) block size from
   % values v.
   assert (numel(blocksize) == 2);

   [i,j] = ndgrid(1:blocksize(1), 1:blocksize(2));

   n     = prod(blocksize);
   N     = numel(v) / n;

   i = repmat(i, [1,N]) + kron(blocksize(1)*(0:N-1), ones(blocksize));
   j = repmat(j, [1,N]) + kron(blocksize(2)*(0:N-1), ones(blocksize));
   S = sparse(i, j, v);
end

%--------------------------------------------------------------------------

function varargout = default_kro(s, krow, krog, Swco, varargin)

   % s = [aqua, liquid, vapour]
   assert (size(s,2) == 3);
   assert (~any(any(s < -sqrt(eps))));

   [krow_s{1 : nargout}] = krow(s(:,2));
   [krog_s{1 : nargout}] = krog(s(:,2));

   sa               = max(s(:,1) - Swco, 0);
   den              = s(:,3) + sa;
   ix               = ~(den > 0);
   varargout{1}     = (s(:,3).*krog_s{1} + sa.*krow_s{1}) ./ den;
   varargout{1}(ix) = 1;

   if nargout > 1,
      dkr = zeros(size(s));

      f     = (krow_s{1} - krog_s{1}) ./ den.^2;
      f(ix) = 0;

      % \partial_{s_a} k_{r,o}
      dkr(:,1)  = s(:,3) .* f;

      % \partial_{s_l} k_{r,o}
      dkr(:,2)  = (s(:,3).*krog_s{2} + sa.*krow_s{2}) ./ den;
      dkr(ix,2) = 0;

      % \partial_{s_v} k_{r,o}
      dkr(:,3)  = - s(:,1) .* f;

      varargout{2} = dkr;
   end
end
