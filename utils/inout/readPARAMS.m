function param = readPARAMS(fn, varargin)
%Read Eclipse or FrontSim data file.
%
% SYNOPSIS:
%   param = readPARAMS(fn)
%
% PARAMETERS:
%   fn - String holding name of (readable) ECLIPSE or FRONTSIM data file.
%        Note: The file is assumed to be a physical(seekable) file on disk,
%        not to be read in through (e.g.) a POSIX pipe.
%
% RETURNS:
%   params -
%
% SEE ALSO:
%   readGRDECL.

%{
#COPYRIGHT#
%}

% $Id: readPARAMS.m 1953 2009-03-31 10:54:12Z bska $

prm  = struct('verbose',  false  ...% emit progress reports
              );
prm = merge_options(prm, varargin{:}); %#ok



[fid, msg] = fopen(fn, 'rt');
if fid < 0, error(msg), end

cartDims    = [-1, -1, -1];
numCell     = -1;
param.TSTEP = [];
param.PSIDE = struct([]);
param.FLUXSIDE = struct([]);
param.NOGRAV   = false;
param.METRIC   = true;
param.WATER    = true;
param.OIL      = true;
param.GAS      = false;
%param.DISGAS   = false;

while ~feof(fid),
  lin = fgetl(fid);
  if lin == -1,
    msg = ferror(fid, 'clear');
    fclose(fid);
    error('readGRDECL:Input:Empty', ...
          'GRDECL file `%s'' is unreadable.\nSystem reports: %s\n', ...
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
    switch kw,
     case {'SPECGRID', 'DIMENS'}
      cartDims = readVector(fid, kw, 3);
      numCell  = prod(cartDims);
      trash    = fgetl(fid);   %#ok
      param.cartDims = cartDims(:).';

     case {'NOGRAV', 'METRIC', 'WATER', 'OIL', 'GAS'}
      param.(kw) = true;

     case 'FIELD'
      param.METRIC = false;

     case {'PSIDE', 'FLUXSIDE'}
      lin    = fgetl(fid);

      tokens = regexp(lin, '\S+', 'match');
      if strcmp(tokens{1}, 'GLOBAL') && ...
         any(strcmp(tokens{2}, ...
                    {'LEFT','RIGHT', ...
                     'FRONT','BACK', ...
                     'TOP','BOTTOM'})),
        vals = str2num(strvcat(tokens{3:end-1})); %#ok

        if numel(vals) == 5,
          param.(kw) = [param.(kw), struct('direction', tokens{2}, ...
                                           'value', vals(1), ...
                                           'range', vals(2:end))];
        else
          error([kw, ':SPEC:Incorrect'], ['Incorrect ', kw, ' Spec found']);
        end
      else
        error([kw, ':SPEC:Incorrect'], ['Incorrect ', kw, ' Spec found']);
      end

     case 'DXV',
      % Error checking?
      checkDim(cartDims, numCell, kw);
      param.(kw) = readVector(fid, kw, cartDims(1));

     case 'DYV',
      checkDim(cartDims, numCell, kw);
      param.(kw) = readVector(fid, kw, cartDims(2));

     case 'DZV',
      checkDim(cartDims, numCell, kw);
      param.(kw) = readVector(fid, kw, cartDims(3));

     case 'DEPTHZ',
      checkDim(cartDims, numCell, kw);
      param.(kw) = readVector(fid, kw, (cartDims(1)+1)*(cartDims(2)+1));

     case 'TSTEP',
      param.TSTEP = [param.TSTEP; readVector(fid, kw, inf)];

     case {'SWAT', 'SGAS', 'PRESSURE', 'RS', 'RV'}
      checkDim(cartDims, numCell, kw);
      param.(kw)  = readVector(fid, kw, numCell);


     otherwise
      %{
      warning('readPARAM:KEYWORD:Unrecognized', ...
              'Input keyword `%s'' is unrecognized\n', kw);
      %}
    end
  end
end
fclose(fid);


%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------

function checkDim(cartDims, numCell, kw)
if numCell < 1 || any(cartDims < 1),
   error('readGRDECL:Input:NoDim', ...
         'GRDECL keyword ''%s'' found before dimension specification', kw);
end
