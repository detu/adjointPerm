function grdecl = readGRDECL(fn, varargin)
%Read Eclipse file assuming corner-point grid format.
%
% SYNOPSIS:
%   grdecl = readGRDECL(fn)
%   grdecl = readGRDECL(fn, 'pn1', pv1, ...)
%
% PARAMETERS:
%   fn - String holding name of (readable) GRDECL specification.
%        Note: The GRDECL specification is assumed to be a physical
%        (seekable) file on disk, not to be read in through (e.g.)
%        a POSIX pipe.
%
%   'pn'/pv - List of 'key'/value pairs designating optional parameters.
%             Currently supported values are
%               - verbose  -- Emit messages to screen while processing.
%                             Default value: FALSE.
%
%               - usecache -- Cache result to mat-file to avoid reading
%                             grid file when it has not changed.
%                             Default value: TRUE.
%
% RETURNS:
%   grdecl - Output structure containing the known, though mostly
%            unprocessed, fields of the GRDECL specification.
%
%            With the exception of keywords regarding wells, stored in
%            'grdecl' structure field 'wells', and 'SPECGRID' whose first
%            three arguments (grid cell dimensions NX, NY, NZ) are stored
%            in 'grdecl' structure field 'cartDims', all currently
%            recognized GRDECL keywords are stored in equally named
%            'grdecl' structure fields. Specifically, GRDECL keyword
%            'COORD' is stored 'grdecl' structure field 'COORD' and so
%            forth.
%
%            The pillar description 'COORD' is stored in an nPillar-by-6
%            array (number of pillars, nPillar == (NX+1)*(NY+1)) of
%            bottom/top coordinate pairs.  Specifically,
%
%              grdecl.COORD(i,1:3) -- (x,y,z)-coordinates of pillar 'i'
%                                     top point.
%              grdecl.COORD(i,4:6) -- (x,y,z)-coordinates of pillar 'i'
%                                     bottom point.
%
%            The currently recognized GRDECL keywords are 'SPECGRID',
%            'COORD', 'ZCORN', 'PORO', 'PERMX', 'PERMY', 'PERMZ', 'ACTNUM',
%            'WELSPECS', 'COMPDAT', 'WCONINJE', 'WCONPROD'.
%
% SEE ALSO:
%   processGRDECL, readWellKW.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-05 14:40:37 +0200 (ma, 05 okt 2009) $
% $Revision: 2948 $


opt = struct('verbose' , false, ...
             'usecache', true , ...
             'cartDims', []   , ...
             'grdecl'  , []);
opt = merge_options(opt, varargin{:});
if isempty(opt.grdecl),
   grdecl = struct;
else
   grdecl = opt.grdecl;
end

[fid, msg] = fopen(fn, 'rt');
if fid < 0, error(msg), end

cartDims = opt.cartDims;
numCell  = [-1, prod(cartDims)];
numCell  = numCell(1 + double(~isempty(cartDims)));

[Vectors, Flags, Wells] = initScheduleVars();

while ~feof(fid),
   lin = fgetl(fid);
   if lin == -1,
      msg = ferror(fid, 'clear');
      fclose(fid);
      error('readGRDECL:Input:Empty', ...
            'GRDECL file ''%s'' is unreadable.\nSystem reports: %s\n', ...
            fn, msg);
   end

   % Loop until next keyword
   kw = regexp(lin, '^[A-Z]+', 'match', 'once');
   while isempty(kw) && ~feof(fid),
      lin = fgetl(fid);
      if lin ~= -1,
         kw = regexp(lin, '^[A-Z]+', 'match', 'once');
      end
   end

   if ~feof(fid),
      dispif(opt.verbose, 'Reading keyword ''%s''\n', kw);
      switch kw,
         case {'SPECGRID', 'DIMENS'},
            t = fscanf(fid, '%f', 3) .';
            if isempty(cartDims),
               cartDims = t;
            elseif ~all(cartDims == t)
               error(msgid('Grid:Initialized'), ...
                     'Attempting to change already defined grid size.');
            end
            numCell  = prod(cartDims);
            trash    = fgetl(fid); %#ok
            grdecl.cartDims = reshape(cartDims, 1, []);

         case 'INCLUDE',
            % checkDim(cartDims, numCell, kw, fid);
            inc_fn_tmp = fscanf(fid, '%s', 1);
            inc_fn = regexp(inc_fn_tmp, '[''"]?([-./\w]+)[''"]?', ...
                            'tokens', 'once');
            inc_fn = inc_fn{1};

            terminated = inc_fn_tmp(end) == '/';
            if inc_fn(end) == '/', inc_fn = inc_fn(1:end-1); end

            % Gobble up keyword-closing '/' character if not already read.
            if ~terminated,
               slash = fscanf(fid, '%s', 1);
               if ~strcmp(slash, '/'),
                  error(msgid('Include:WrongfulTermination'), ...
                        'INCLUDE keyword not correctly terminated.');
               end
            end

            inc_fn(inc_fn == '/') = filesep;
            if inc_fn(1) ~= filesep,
               % Translate relative pathname to absolute pathname.
               inc_fn = fullfile(fileparts(fn), inc_fn);
            end

            dispif(opt.verbose, ' -> ''%s''\n', inc_fn);
            grdecl = readGRDECL(inc_fn, 'cartDims', cartDims, ...
                                    'verbose', opt.verbose,       ...
                                    'usecache', opt.usecache,     ...
                                    'grdecl',grdecl);
            dispif(opt.verbose, ' <- ''%s''\n', inc_fn);

         case 'COORD',
            checkDim(cartDims, numCell, kw, fid);
            grdecl.COORD = readVector(fid, kw, ...
                                      6 * prod(cartDims(1:2) + 1));

         case 'ZCORN',
            checkDim(cartDims, numCell, kw, fid);
            grdecl.ZCORN = readVector(fid, kw, 8 * numCell);

         case {'PORO',                       ...
               'PERMX',  'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY',  'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ',  ...
               'ACTNUM', 'SATNUM', 'PVTNUM', ...
               'MULTX',  'MULTX-',           ...
               'MULTY',  'MULTY-',           ...
               'MULTZ',  'MULTZ-',           ...
               'EQLNUM'
               %'FLUXNUM,
               },
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(regexprep(kw, '\W', '_')) = readVector(fid, kw, numCell);

         case {'COORDX', 'COORDY', 'COORDZ'},
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(kw) = readVector(fid, kw, prod(cartDims+1));

         case {'NOGRAV', 'METRIC', 'WATER', 'OIL', 'GAS'},
            if ~Flags.(kw).set, Flags.(kw).set = true; end
            Flags.(kw).val = true;

         case 'FIELD',
            if ~Flags.METRIC.set, Flags.METRIC.set = true; end
            Flags.METRIC.val = false;

         case {'PSIDE', 'FLUXSIDE'},
            lin    = fgetl(fid);
            tokens = regexp(lin, '\S+', 'match');

            if strcmp(tokens{1}, 'GLOBAL') && ...
               any(strcmp(tokens{2}, ...
                          {'LEFT', 'RIGHT', ...
                           'FRONT', 'BACK', ...
                           'TOP', 'BOTTOM'})),
               vals = str2num(strvcat(tokens{3:end-1})); %#ok

               if numel(vals) == 5,
                  if ~Vectors.(kw).set, Vectors.(kw).set = true; end
                  Vectors.(kw).val = [Vectors.(kw).val, ...
                                      struct('direction', tokens{2}, ...
                                             'value', vals(1),       ...
                                             'range', vals(2:end))];
               else
                  error([kw, ':SPEC:Incorrect'], ['Incorrect ', kw, ' Spec found']);
               end
            else
               error([kw, ':SPEC:Incorrect'], ['Incorrect ', kw, ' Spec found']);
            end

         case {'DXV', 'DYV', 'DZV'},
            % Error checking?
            checkDim(cartDims, numCell, kw, fid);
            ix          = strcmp(kw, {'DXV', 'DYV', 'DZV'});
            grdecl.(kw) = readVector(fid, kw, cartDims(ix));

         case 'DEPTHZ',
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(kw) = readVector(fid, kw, prod(cartDims(1:2) + 1));

         case 'TSTEP',
            if ~Vectors.TSTEP.set, Vectors.TSTEP.set = true; end
            Vectors.TSTEP.val = [Vectors.TSTEP.val; readVector(fid, kw, inf)];

         case {'SWAT', 'SGAS', 'PRESSURE', 'RS', 'RV'},
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(kw)  = readVector(fid, kw, numCell);

         case {'WELSPECS', 'COMPDAT', 'WCONINJE', 'WCONPROD', 'WCONHIST'},
            if ~Wells.set, Wells.set = true; end
            Wells.val = readWellKW(fid, Wells.val, kw);

         case {'EQUALS', 'MULTIPLY'},
            grdecl = readOperatorModField(fid,grdecl,kw);

         case {'COPY'},
            lin = fgetl(fid);
            count = 1;
            while isempty(regexp(lin, '^/', 'once')),
               % Skip blank lines and comments.
               if(~(isempty(lin) || ~isempty(regexp(lin, '^--', 'match'))))
                  split = regexp(lin, '(\w\.*\-*\**)+', 'match');
                  if(~length(split)==2)
                     error('Wrong line in COPY')
                  end
                  name1 = char(split(1));
                  name2 = char(split(2));
                  if(isfield(name2,grdecl)),
                     error('Error: try to copy to exising field')
                  else
                     grdecl.(name2) = grdecl.(name1);
                  end
               end
               lin = fgetl(fid);
               count =count +1;
            end

         case {'FAULTS'},
            grdecl = readFaults(fid,grdecl);

         case {'MULTFLT'},
            grdecl = readMultfelt(fid,grdecl);

         case {'EQUIL'},
            grdecl = readEQUIL(fid,grdecl);

         case {'NOECHO', 'PROPS', 'MESSAGES', 'START', 'EQLDIMS', ...
               'REGDIMS', 'WELLDIMS', 'TABDIMS', 'VFPIDIMS',      ...
               'VFPPDIMS', 'FAULTDIM', 'PIMTDIMS', 'NSTACK',      ...
               'UNIFIN', 'UNIFOUT', 'OPTIONS', 'NEWTRAN',         ...
               'SCHEDULE','SOLUTION'},
            dispif(opt.verbose > 1, 'Keyword ''%s'' ignored.\n', kw);

         otherwise,
            dispif(opt.verbose > 1, ...
                   'Input keyword ''%s'' is unrecognized\n', kw);
      end
   end
end
fclose(fid);
if isfield(grdecl, 'ACTNUM'),
   grdecl.ACTNUM = int32(grdecl.ACTNUM);
end
grdecl = mergeScheduleVars(grdecl, Vectors, Flags, Wells);

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function well = initDefaultWell()
% Initialize the default well structure
% NOTE:
% * all fields ending with '_vec' are vectors with values for each COMPDAT.
% * for comments define: rate_trg = rate target, up_lim = upper limit.

   well = struct( ...
         'name',       '', ...
         'I',          0,  ...
         'J',          0,  ...
         'I_vec',      0, ...
         'J_vec',      0, ...
         'K_vec_u',    -1, ... %Location of upper connecting block.
         'K_vec_l',    -1, ... %Location of lower connecting block.
         'state_vec',   1, ... %OPEN/SHUT flag of connection, 1 = 'OPEN'.
         'diam_vec',   0.25, ...
         'transm_vec', -1, ... %Transm. factor of connection (WI).
         'Kh_vec',     -1, ... %Permeability thickness.
         'skin_vec',   0.0,... %Skin factor of connection.
         'dir_vec',    'Z',... %Dir of wellblock penetration: X, Y or Z.
         'state',      'OPEN', ... %OPEN or SHUT well.
         'control',    '', ... %rat_trg type or BHP.
         'oil_rate',   [], ... %rat_trg/up_lim (prod).
         'water_rate', [], ... %rat_trg/up_lim (prod).
         'gas_rate',   [], ... %rat_trg/up_lim (prod).
         'liquid_rate',[], ... %rat_trg/up_lim (prod) or surface rate (inj).
         'vol_rate',   [], ... %rat_trg/up_lim (reservoir fluid vol-rate).
         'bhp',        [], ... %bhp target or lower limit.
         'type',       '', ... %Injection/Production.
         'inj_type',   '', ... %WAT, GAS, OIL.
         'compdat_count', 0);  %Count number of compdats.

%--------------------------------------------------------------------------

function checkDim(cartDims, numCell, kw, fid)
if isempty(cartDims) || numCell < 1 || any(cartDims < 1),
   % Don't leave open fd's in MATLAB's workspace when erroring out.
   fclose(fid);
   error('readGRDECL:Input:NoDim', ...
         'GRDECL keyword ''%s'' found before dimension specification', kw);
end

%--------------------------------------------------------------------------

function [Vectors, Flags, Wells] = initScheduleVars()
[Vectors.PSIDE,    ...
 Vectors.FLUXSIDE, ...
 Vectors.TSTEP]        = deal(struct('set', false, 'val', []));

[Flags.METRIC, ...
 Flags.NOGRAV, ...
 Flags.WATER,  ...
 Flags.OIL,    ...
 Flags.GAS]            = deal(struct('set', false, 'val', []));

Wells                  = struct('set', false, ...
                                'val', initDefaultWell());

%--------------------------------------------------------------------------

function grdecl = mergeScheduleVars(grdecl, Vectors, Flags, Wells)
if isempty(grdecl), grdecl = struct(); end
fn1 = @(s,n) reshape(n(cellfun(@(f) s.(f).set, n)), 1, []);
fn  = @(s  ) fn1(s, fieldnames(s));

for nm = fn(Vectors), nm1 = nm{1}; grdecl.(nm1) = Vectors.(nm1).val; end
for nm = fn(Flags  ), nm1 = nm{1}; grdecl.(nm1) = Flags  .(nm1).val; end

if Wells.set, grdecl.wells = Wells.val(2:end); end

%--------------------------------------------------------------------------

%{
function grdecl = mergeIncludeFile(grdecl, inc_grdecl)
if isempty(grdecl), grdecl = struct(); end
for nm = reshape(fieldnames(inc_grdecl), 1, []),
   nm1 = nm{1};
   if any(strcmp(nm1, {'FLUXSIDE', 'PSIDE', 'TSTEP'})),
      grdecl.(nm1) = [grdecl.(nm1), inc_grdecl.(nm1)];
   elseif any(strcmp(nm1,{'MULTIPLY', 'EQUALS'})),
       for nmm = reshape(fieldnames(inc_grdecl.(nm1)), 1, []),
           nm2 = nmm{1};
           if(~isfield(nm2,grdecl.(nm1))),
               grdecl.(nm1).(nm2) = inc_grdecl.(nm1).(nm2);
           else
               grdecl.(nm1).(nm2).value =[grdecl.(nm1).(nm2).value;inc_grdecl.(nm1).value];
               grdecl.(nm1).(nm2).region =[grdecl.(nm1).(nm2).region;inc_grdecl.(nm1).region];
           end
       end
   else
       grdecl.(nm1) = inc_grdecl.(nm1);
   end
end
%}
