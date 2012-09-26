function v = readVector(fid, field, nel)
%Input vector of floating point numbers from GRDECL file.
%
% SYNOPSIS:
%   v = readVector(fid, field, nel)
%
% PARAMETERS:
%   fid   - File identifier (as defined by FOPEN) of GRDECL file open for
%           reading.  Assumed to point to a seekable (i.e. physical) file
%           on disk, not (e.g.) a POSIX pipe.
%
%   field - Name (string) identifying the GRDECL keyword (field) currently
%           being processed.  Used for error identification/messages only.
%
%   nel   - Number of elements to read from input stream.  As a special
%           case, the caller may pass nel==INF (or nel=='inf') to read as
%           much as possible.  In this case it is imperative that the input
%           vector be terminated by a '/' character.
%
% RETURNS:
%   v     - Vector, length 'nel', of floating point numbers defining the
%           contents of GRDECL keyword 'field'.  If nel==INF, NUMEL(v) is
%           the number of vector elements read before the terminating slash
%           character.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-18 10:38:11 +0200 (fr, 18 sep 2009) $
% $Revision: 2770 $

% Simple vector read.  Assume the vector is given explictly by listing all
% of its values separated by whitespace, possibly interspersed by comments.
%

if ischar(nel) && strcmpi(nel, 'inf'), nel = inf; end
scanner = getScanner();
opt = getScanOpt(nel);
cs  = opt{end};   % CommentStyle -- relies on knowledge of 'getScanOpt'.
C   = scanner(fid, '%f', opt{:});

pos = ftell(fid);
lin = fgetl(fid);
while ischar(lin) && (isempty(lin) || matches(lin, ['^\*|^\s+$|^', cs])),
   if strmatch('*', lin),
      % Repeat description of the form
      %      N*value
      % with 'N' a positive integer.  There must be no whitespace
      % immediately before or after the literal asterisk ('*').
      %
      % This format allows considerable economy of file space, e.g. when
      % listing uniform porosity values such as
      %
      %   PORO
      %     5000*0.3 /
      %
      % for a vector of five thousand repetitions of the value 0.3.  The
      % MATLAB equivalent is
      %
      %   PORO = REPMAT(0.3, [5000, 1]);
      %
      % This format is most likely encountered when processing GRDECL
      % keywords such as 'DXV', 'PERMX', 'PORO', 'ACTNUM' and, possibly,
      % 'ZCORN' and 'TSTEP'.
      %

      % Validate repeat description format
      fseek(fid, pos - 1, 'bof');
      repdesc = fscanf(fid, '%c', 3);
      if numel(repdesc) ~= 3 || ~matches(repdesc, '\d\*\d'),
         error('readGRDECL:RepeatDescr:Malformed',                    ...
               'Incorrect repeat description detected in keyword %s', ...
               field);
      end

      % Validate repeat count value.
      nrep = C{1}(end);
      if ~isnumeric(nrep)             || ...
         nrep < 1                     || ...
         nrep > nel - numel(C{1}) + 1 || ...
         nrep ~= fix(nrep),
         error('readGRDECL:RepeatCount:OutOfBounds',             ...
               ['Unexpected repeat count in keyword %s.\n',      ...
                'Expected integral repeat count in [1 .. %d], ', ...
                'but found %f.'], field, nel - numel(C{1}) + 1, nrep);
      end

      % Append repeated value to vector.
      [val, count, errmsg, nextix] = sscanf(lin(2:end), '%f', 1);
      C   = {[C{1}(1:end-1); repmat(val, [nrep, 1])]};
      nv  = numel(C{1});

      % Reposition reader at first whitespace or '/' character whereupon
      % regular vector reading may resume.
      fseek(fid, pos + nextix, 'bof');

      opt = getScanOpt(nel, nv);
      C1 = scanner(fid, '%f', opt{:});
      C  = {[C{1}(:); C1{1}(:)]};
   end
   pos = ftell(fid);
   lin = fgetl(fid);
end

iseof           = feof(fid);
[errmsg, errno] = ferror(fid);
readAll = (~isinf(nel) && (numel(C{1}) == nel)) | isinf(nel);
if ischar(lin) && matches(lin, '^\s*/') && readAll,
   v = C{1};
elseif ischar(lin) && matches(lin, '^\s*/'),
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':PrematureTermination'],    ...
         ['Detected termination character ''/'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.'], pos, field);
elseif ischar(lin) && matches(lin, '^[A-Z]+'),
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':OtherKeyword'],           ...
         ['Encountered non-numeric data ''%s'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.\n',           ...
          'Missing slash (/) in input?'],                          ...
          lin, pos, field);
elseif ischar(lin) && matches(lin, '^[a-z!;=]+'),
   error(['readGRDECL:Vector_', field, ':Unexpected'],             ...
         ['Encountered non-numeric data ''%s'' at position %d.\n', ...
          'Scanning of keyword ''%s'' not completed.'],            ...
          lin, pos, field);
elseif iseof,
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':EOF'],               ...
         'End of file before complete read of keyword ''%s''.', field);
elseif errno ~= 0,
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':SystemError'], ...
         ['Read error while processing keyword %s.\n',  ...
          'Input system reports ''%s''.'], field, errmsg);
else
   fclose(fid);
   error(['readGRDECL:Vector_', field, ':Malformed'],     ...
         ['Unexpected input while reading keyword %s.\n', ...
          'Found ''%s'' at position %d.'], field, lin, pos);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function scan = getScanner()
if exist('textscan', 'builtin'),
   scan = @textscan;
else
   scan = @myscanner;
end

%--------------------------------------------------------------------------

function opt = getScanOpt(nel, varargin)
opt = {'CommentStyle', '--'};
if nargin > 1 && isnumeric(varargin{1}),
   nel = nel - varargin{1};
end
if ~isinf(nel),
   opt = {nel, opt{:}};
end

%--------------------------------------------------------------------------

function b = matches(str, pat)
b = ~isempty(regexp(str, pat, 'once'));

%--------------------------------------------------------------------------

function C = myscanner(fid, format, varargin)
%Scanner function supporting a subset of TEXTSCAN features.

% This function is needed when TEXTSCAN is unavailable (e.g., Octave or old
% (pre-R14) releases of M).
%
% Algorithm:
%   1) Attempt to read as much data as possible using FSCANF.  This handles
%      most common cases (all data listed separately, no comments).
%   2) If we encounter a comment, as identified by 'opt.CommentStyle', then
%      scan individual lines until all comments have been scanned.
%   3) If we encounter 'end-of-record' (i.e., a '/' character), then
%      terminate all scanning and trust caller to detect any errors.

% Setup:
%   Define comment style (usually '--' to exclude a single line of input)
%   Define termination strings (regexp's).
nel = inf;
if mod(numel(varargin), 2) == 1,
   if readnum_ok(varargin{1}),
      nel = varargin{1};
   end
   varargin = { varargin{2 : end} };
end
opt = struct('CommentStyle', '--');
opt = merge_options(opt, varargin{:});
if iscell(opt.CommentStyle),
   assert (numel(opt.CommentStyle) == 1);
   opt.CommentStyle = opt.CommentStyle{1};
end
assert (ischar(opt.CommentStyle));
patt_comment   = ['^', opt.CommentStyle];
patt_terminate = '^\*|^/';

% Algorithm step 1.  Consume as much data as possible.
[v, cnt] = fscanf(fid, format, nel);
n = 0;
C = reshape(v, [], 1);
read_complete = n + cnt == nel;
while ~read_complete,
   % Short read.  Encountered comment or "premature" termination (usually
   % when ISINF(nel)).  Must either consume a (block of) comment line(s)
   % and all subsequent data or detect end-of-record in which case we need
   % to terminate.

   % Algorithm:
   %   1) Record position in file to allow caller to inspect input stream
   %      if we need to terminate.
   %   2) Inspect the next few character to determine if this is a comment
   %      or 'end-of-record'.
   %   3) Act accordingly.
   pos = ftell(fid);
   lin = fgetl(fid);
   if strmatch(opt.CommentStyle, '--') && matches(lin, '^-'),
      % CommentStyle is '--' (i.e., following Eclipse convention).  FSCANF
      % consumed the first '-' character thinking that the next token would
      % be (negative) floating point number and failed when the second '-'
      % character was encountered.  Back up one character to allow comment
      % matcher to determine that this line is indeed a comment.
      fseek(fid, pos-1, 'bof');
      pos = ftell(fid);
      lin = fgetl(fid);
   end
   n = n + cnt;
   if matches(lin, patt_comment),
      % This is a comment.  The above FGETL has already consumed all
      % character up to (and including) EOL.  Attempt to read new data with
      % regular FSCANF approach.  This will fail if the next line is
      % another comment, but then the enclosing WHILE will bring us back
      % here on the next iteration...
      [v, cnt] = fscanf(fid, format, nel - n);
      C = [C; reshape(v, [], 1)];  %#ok
      read_complete = n + cnt == nel;
   elseif matches(lin, patt_terminate),
      % This is a termination character (i.e., the '/' character) or a
      % repeat character (i.e., the '*' character).  Back up to beginning
      % of 'lin' to allow caller to rediscover the special character and
      % take appropriate action.  Our work here is done.
      fseek(fid, pos, 'bof');
      read_complete = true;
   else
      fclose(fid);
      error('readVector:Input:Unexpected', ...
           ['Unexpected input ''%s'' encountered at position ''%d''', ...
            ' whilst scanning file.'], lin, pos);
   end
end

% Convert to CELL to conform to semantics of TEXTSCAN.
C = { C };

%--------------------------------------------------------------------------

function b = readnum_ok(n)
b = isnumeric(n) && ~isempty(n) && isfinite(n);
