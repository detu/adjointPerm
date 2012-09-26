function [pvt, density, miscible, info] = initBlackOilPVT(T, varargin)
%Construct Black Oil PVT model evaluation function.
%
% SYNOPSIS:
%   [pvt, density, miscible, info] = initBlackoilPVT(T)
%   [pvt, density, miscible, info] = initBlackoilPVT(T, 'pn1', pv1, ...)
%
% PARAMETERS:
%   T       - ECLIPSE PVT table structure as defined by (e.g.) function
%             'readpvt'.  It must contain valid fields '.pvto(.pvdo)',
%             '.pvtg(.pvdg)', '.pvtw', and '.density'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Verbose -- Whether or not to emit informational messages
%                            while interpolating PVT tables.
%                            Logical.  Default value = FALSE.
%
% RETURNS:
%   pvt      - Function which evaluates PVT data for the specific Black Oil
%              model.  Specifically, given n-by-1 pressures p and n-by-np
%              surface volumes z, the call
%
%                  [c, rho, mu, u, R] = pvt(p, z)
%
%              computes, respectively, n-by-np phase compressibilities (c),
%              densities (rho), viscosities (mu), surface volumes (u), and
%              mixing ratio (R).
%
%   density  - An np vector of surface densities of each phase/component.
%
%   miscible - An np boolean vector indicating whether or not the
%              respective phase/component mixes in the other phases.
%
%   info     - A textual representation of the specific fluid properties
%              (i.e., surface densities and miscibility).
%
% EXAMPLE:
%   tab = readpvt(fullfile(ROOTDIR, 'params', 'fluid', ...
%                          'Data', 'miscibleoil.txt'));
%   [pvt, density, miscible, info] = initBlackOilPVT(tab, 'Verbose', true);
%
% SEE ALSO:
%   readpvt, initBlackoilFluid, initBlackOilRelPerm.

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
      if sum(illegal)==1, error('Unknown phase specified');
      else                error('Unknown phases specified'); end
   end

   if numel([A,L,V]) == 0, error('Huh!?, No phases?'); end

   % Get number of pvt regions
   num_pvt_regions = 1;
   if isfield(T, 'pvtnum'), num_pvt_regions = max(T.pvtnum);   end
   if num_pvt_regions > 1, error ('PVTNUM not supported yet'); end

   check = @(T, f) isfield(T,f) && ...
      ~(iscell(T.(f)) && numel(T.(f))~=num_pvt_regions);


   % Surface densities
   if check(T, 'density'), d = T.density(1,:);
      assert(numel(d)==3);
      if ~isempty(A),  density(A) = d(2);      end
      if ~isempty(L),  density(L) = d(1);      end
      if ~isempty(V),  density(V) = d(3);      end
   else
      error (['PVT table must contain keyword "density" with one table ',...
         'per pvt region']);
   end

   miscible =false(1,max([A,L,V]));

   % Water phase pvt function.
   if ~isempty(A),
      if check(T, 'pvtw'),
         pvtfun{A} = @(p,z) incompressiblewater(T.pvtw(1,:), p, z(:,A), z(:,L), z(:,V), opt.verbose);
      else
         error(id('WaterPhaseData:NotPresent'), ...
            'No data or too much data for water phase.\n');
      end
   end

   % Gas phase pvt function
   if ~isempty(V),
      if ~xor(check(T, 'pvtg'), check(T, 'pvdg'))
         error(id('GasPhaseData:DataError'), ...
            ['No data or too much data for gas phase. One and only one of the\n',...
            'fields "pvtg" and "pvdg" must be present and the number of tables\n',...
            'must match the number of pvt regions']);

      end

      if check(T,'pvtg'),
         if isempty(L),
            error('Miscible vapor phase with no oil phase is impossible.');
         end
         pvtfun{V}=@(p,z) wet_gas(T.pvtg, p, z(:,A), z(:,L), z(:,V), opt.verbose);
         miscible(V) = true;

      elseif check(T,'pvdg'),
         pvtfun{V}=@(p,z) dry_gas(T.pvdg, p, z(:,A), z(:,L), z(:,V), opt.verbose);
      end
   end

   % Oil phase pvt function
   if ~isempty(L),
      if ~xor(check(T, 'pvto'), check(T, 'pvdo'))
         error(id('OilPhaseData:DataError'), ...
            ['No data or too much data for oil phase. One and only one of the\n',...
            'fields "pvto" and "pvdo" must be present and the number of tables\n',...
            'must match the number of pvt regions']);
      end

      if check(T,'pvto'),
         if isempty(V),
            error('Miscible liquid phase with no gas phase is impossible.');
         end
         pvtfun{L}=@(p,z) live_oil(T.pvto, p, z(:,A), z(:,L), z(:,V), opt.verbose);
         miscible(L) = true;

      elseif check(T,'pvdo'),
         pvtfun{L}=@(p,z) dead_oil(T.pvdo, p, z(:,A), z(:,L), z(:,V), opt.verbose);
      end
   end

   if opt.verbose, disp(getInfo(density, miscible, A, L, V)); end

   % ==========================================================/
   %
   % This function is the return value from initBlackOilPVT.
   %
   % For convenience, this function binds the PVT tables, and takes only
   % pressure (p) and surface volumes (z) as input arguments.
   %
   % Alternativley, the function takes a single string argument to lookup
   % info, phase surfacedensities, miscibility and names.
   %
   function [c, rho, mu, u, R] = blackoilPVT(p, z)
      %[C, RHO, MU, U, FRAC] = blackoilPVT(varargin)

      if numel(p) ~= size(z, 1),
         error('There must be one pressure for each row in mass table');
      end
      if size(z, 2) ~= numel([A, L, V]),
         error ('There must be one column in mass table for each phase');
      end
      % Convert units of from SI to FIELD or METRIC, i.e.,
      % convert p from Pascal to Psia or Bar.
      if isfield(T, 'field'),
         p = convertTo(p, psia());
      else
         p = convertTo(p, barsa());
      end

      [B, dB, R, dR, mu, P] = deal(zeros([numel(p), numel(pvtfun)]));
      for i = 1 : numel(pvtfun),
         [B(:,i), dB(:,i), ...
          R(:,i), dR(:,i), mu(:,i), P(:,i)] = pvtfun{i}(p, z);
      end

      % If the vapor phase is undersaturated, the liquid phase is missing.
      % Consequently, no gas can be dissolved in liquid.
      R(~P(:,V) | ~(z(:,L)>0), L) = 0;

      % Likewise, if the liquid phase is undersaturated, the vapor phase is
      % missing.  Consequently, no oil is allowed to evaporate.
      R(~P(:,L) | ~(z(:,V)>0), V) = 0;

      [c, rho, u] = convertBlackoil(z, B, dB, R, dR, density, {A,L,V});

      % Convert units from FIELD or METRIC to SI
      if isfield(T, 'field'),
         c   = convertFrom(c,   1/psia());
         rho = convertFrom(rho, 1); %??
         assert (false);
      else
         c   = convertFrom(c,   1/barsa());
      end

      mu = convertFrom(mu, centi*poise());
   end

   pvt = @blackoilPVT;
   info = getInfo(density, miscible, A, L, V);
end

%--------------------------------------------------------------------------

function s = getInfo(density, miscible, A, L, V)
   s = '';
   if A, s = [s,sprintf('Surface density for water: %8.2f\n',   density(A))]; end
   if L, s = [s,sprintf('Surface density for oil:   %8.2f\n',   density(L))]; end
   if V, s = [s,sprintf('Surface density for gas:   %8.2f\n\n', density(V))]; end

   if A,
      if miscible(A),s = [s,sprintf('Aquaic phase:  miscible\n')];
      else           s = [s,sprintf('Aquaic phase:  immiscible\n')]; end
   end

   if L,
      if miscible(L),s = [s,sprintf('Liquid phase:  miscible\n')];
      else           s = [s,sprintf('Liquid phase:  immiscible\n')]; end
   end

   if V,
      if miscible(V),s = [s,sprintf('Vapour phase:  miscible\n')];
      else           s = [s,sprintf('Vapour phase:  immiscible\n')]; end
   end
end

%--------------------------------------------------------------------------

function [B,dB,R,dR,mu,P] = incompressiblewater(pvtw, p, zA, zL, zV, verbose) %#ok
   B   = repmat(1.0,     size(p));
   dB  = repmat(0.0,     size(p));
   R   = repmat(0.0,     size(p));
   dR  = repmat(0.0,     size(p));
   mu  = repmat(pvtw(4), size(p));
   P   = true(size(p));
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = dry_gas(pvdg, p, zA, zL, zV, verbose) %#ok
%----------------------------------------
   [B, dB, mu] = immiscible_table(pvdg, p, verbose);
   R  = zeros([numel(p),1]);
   dR = zeros([numel(p),1]);
   P  = true(size(p));
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = wet_gas(pvtg, p, zA, zL, zV, verbose)%#ok
%----------------------------------------
   [B, dB, R, dR, mu, P] = miscible_gas(pvtg, p, zL ./ zV, verbose);
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = dead_oil(pvdo, p, zA, zL, zV, verbose) %#ok
%----------------------------------------
   [B, dB, mu] = immiscible_table(pvdo, p, verbose);
   R  = zeros([numel(p),1]);
   dR = zeros([numel(p),1]);
   P  = true(size(p));
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = live_oil(pvto, p, zA, zL, zV, verbose) %#ok
%----------------------------------------
   [B, dB, R, dR, mu, P] = miscible_oil(pvto, p, zV ./ zL, verbose);
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, mu] = immiscible_table(tab, p, verbose) %#ok
%----------------------------------------
   B  = interpTable (tab(:,1), tab(:,2), p);
   dB = dinterpTable(tab(:,1), tab(:,2), p);
   mu = interpTable (tab(:,1), tab(:,3), p);
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = miscible_oil(tab, p, volratio, verbose)%#ok
%----------------------------------------
   % Declare size
   sz = size(p);
   B  = nan(sz);
   dB = nan(sz);
   dR = nan(sz);
   mu = nan(sz);
   P  = true(sz);
   % Find row numbers such that each table section i
   % spans rows row(i):row(i+1)-1
   k       = find(~isnan(tab(:,1)));
   row     = [find(~isnan(tab(:,1))); length(tab(:,1))+1];

   % --- saturated case
   R       = interpTable (tab(k,2), tab(k,1), p);
   maxR    = volratio;
   %  R       = min (R, maxR);

   % If R < maxR, we have saturated case
   is      = find(R <  maxR);

   mu(is) = interpTable (tab(k,2), tab(k,4), p(is));
   dR(is) = dinterpTable(tab(k,2), tab(k,1), p(is));
   B (is) = interpTable (tab(k,2), tab(k,3), p(is));
   dB(is) = dinterpTable(tab(k,2), tab(k,3), p(is));

   % undersaturated case
   iu      = find(R >= maxR);
   if any(iu)

      R(iu)  = maxR(iu);
      dR(iu) = 0;
      %
      r      = maxR(iu);
      p      = p(iu);

      % group points together that refer to same table section
      delim   = [0;tab(k,1);inf];
      nbin    = numel(delim)-1;
      [n,bin] = histc  (r, delim);
      groups  = unique (bin);

      % Compute Weighted average of Bo, dBo and mu for each rl
      for i=1:length(groups)
         g      = groups(i);
         I      = bin==g;
         p1     = p(I);
         iuI    = iu(I);


         if g==1
            % Extrapolate from first table section
            dispif(false, 'Extrapolate from first table section\n');
            i2      = row(g):row(g+1)-1;
            mu(iuI) = interpTable (tab(i2,2), tab(i2,4), p1);
            B (iuI) = interpTable (tab(i2,2), tab(i2,3), p1);
            dB(iuI) = dinterpTable(tab(i2,2), tab(i2,3), p1);

         elseif g==nbin
            % Extrapolate from last table section
            dispif(false, 'Extrapolate from last table section\n');
            i1      = row(g-1):row(g)-1;
            mu(iuI) = interpTable (tab(i1,2), tab(i1,4), p1);
            B (iuI) = interpTable (tab(i1,2), tab(i1,3), p1);
            dB(iuI) = dinterpTable(tab(i1,2), tab(i1,3), p1);

         else

            % Interpolate between table sections
            tab_r   = [tab(row(g-1),1), tab(row(g),1)];
            w       = [tab_r(2)-r(I), r(I)-tab_r(1)]./diff(tab_r);

            i1      = row(g-1):row(g)-1;
            mu1     = interpTable (tab(i1,2), tab(i1,4), p1);
            B1      = interpTable (tab(i1,2), tab(i1,3), p1);
            dB1     = dinterpTable(tab(i1,2), tab(i1,3), p1);

            i2      = row(g):row(g+1)-1;
            mu2     = interpTable (tab(i2,2), tab(i2,4), p1);
            B2      = interpTable (tab(i2,2), tab(i2,3), p1);
            dB2     = dinterpTable(tab(i2,2), tab(i2,3), p1);

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = w(:,1) .* B1  + w(:,2) .* B2;
            dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;
         end
      end
   end
   P(iu) = false;
end

%--------------------------------------------------------------------------

%----------------------------------------
function [B, dB, R, dR, mu, P] = miscible_gas(tab, p, volratio, verbose)%#ok
%----------------------------------------
   % Declare size
   sz  = size(p);
   B  = nan(sz);
   dB = nan(sz);
   dR = nan(sz);
   mu = nan(sz);
   P  = true(sz);

   % Find row numbers such that each table section i
   % spans rows row(i):row(i+1)-1
   k       = find(~isnan(tab(:,1)));
   row     = [find(~isnan(tab(:,1))); length(tab(:,1))+1];

   % --- saturated case
   R       = interpTable (tab(k,1), tab(k,2), p);
   maxR    = volratio;
   %R      = min (R, maxR);

   % If R < maxR, we have saturated case
   is      = find(R <  maxR);

   mu(is) = interpTable (tab(k,1), tab(k,4), p(is));
   dR(is) = dinterpTable(tab(k,1), tab(k,2), p(is));
   B (is) = interpTable (tab(k,1), tab(k,3), p(is));
   dB(is) = dinterpTable(tab(k,1), tab(k,3), p(is));

   % undersaturated case
   iu      = find(R >= maxR);
   if any(iu)

      R(iu)  = maxR(iu);
      dR(iu) = 0;
      %
      r      = R(iu); %maxR(iu);
      p      = p(iu);

      % group points together that refer to same table section
      delim   = [0;tab(k,1);inf];
      nbin    = numel(delim)-1;
      [n,bin] = histc  (p, delim);
      groups  = unique (bin);

      % Compute Weighted average of Bo, dBo and mu for each rl
      for i=1:length(groups)
         g      = groups(i);
         I      = bin==g;
         r1     = r(I);
         p1     = p(I);
         iuI    = iu(I);


         if g==1
            % Extrapolate from first table section
            dispif(false, 'Extrapolate from first table section\n');
            i2      = row(g):row(g+1)-1;
            mu(iuI) = interpTable (tab(i2,2), tab(i2,4), r1);
            B (iuI) = interpTable (tab(i2,2), tab(i2,3), r1);
            %dB(iuI) = dinterpTable(tab(i2,2), tab(i2,3), r1); % GALT

            i3      = row(g+1):row(g+2)-1;
            B3      = interpTable (tab(i3,2), tab(i3,3), r1);
            dB(iuI) = (B3-B(iuI))./(tab(row(g+1),1)- tab(row(g),1));

         elseif g==nbin
            % Extrapolate from last table section
            dispif(false, 'Extrapolate from last table section\n');
            i1      = row(g-1):row(g)-1;
            mu(iuI) = interpTable (tab(i1,2), tab(i1,4), r1);
            B (iuI) = interpTable (tab(i1,2), tab(i1,3), r1);
            %dB(iuI) = dinterpTable(tab(i1,2), tab(i1,3), r1); % GALT

            i0      = row(g-2):row(g-1)-1;
            B0      = interpTable (tab(i0,2), tab(i0,3), r1);
            dB(iuI) = (B(iuI)-B0)./(tab(row(g-1),1)- tab(row(g-2),1));
         else

            % Interpolate between table sections
            tab_p   = [tab(row(g-1),1), tab(row(g),1)];
            w       = [tab_p(2)-p1, p1-tab_p(1)]./diff(tab_p);

            i1      = row(g-1):row(g)-1;
            mu1     = interpTable (tab(i1,2), tab(i1,4), r1);
            B1      = interpTable (tab(i1,2), tab(i1,3), r1);
            %dB1     = dinterpTable(tab(i1,2), tab(i1,3), r1); % GALT

            i2      = row(g):row(g+1)-1;
            mu2     = interpTable (tab(i2,2), tab(i2,4), r1);
            B2      = interpTable (tab(i2,2), tab(i2,3), r1);
            %dB2     = dinterpTable(tab(i2,2), tab(i2,3), r1); % GALT

            mu(iuI) = w(:,1) .* mu1 + w(:,2) .* mu2;
            B (iuI) = w(:,1) .* B1  + w(:,2) .* B2;
            %dB(iuI) = w(:,1) .* dB1 + w(:,2) .* dB2;

            % Use simple FD approximation
            dB(iuI) = (B2-B1)./(tab_p(2)-tab_p(1));
         end
      end
   end
   P(iu) = false;
end

%
% The PVTG tables are slightly different from the PVTO tables:
% Firstly, the Po(Pg) and Rl (Rv) values are in column 2(1) and
% 1 (2).  Secondly, while Rl increase with increasing oil pressure,
% Rv decrease with increasing gas pressure.

%--------------------------------------------------------------------------

function s = id(s)
   s = ['initBlackOilPVT:', s];
end
