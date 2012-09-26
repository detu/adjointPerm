function Table = readpvt(filename, varargin)
%Read Black Oil PVT data from ECLIPSE parameter file into table.
%
% SYNOPSIS:
%   Table = readpvt(filename)
%   Table = readpvt(filename, 'pn', pv, ...)
%
% PARAMETERS:
%   filename - Name of file from which to read PVT data.
%
%   'pn'/pv  - List of 'key'/value paris defining optional parameters.  The
%              supported options are:
%                - verbose -- Whether or not to emit PVT table section
%                             names as new sections are encountered.
%                             Logical.  Default value = FALSE.
%
% RETURNS:
%   Table    - Struct of Black Oil Tables.  Currently, the following ECLIPSE
%              keywords are supported:
%
%              PVTO    - properties of live oil (with dissolved gas)
%              PVDO    - properties of dead oil (no dissolved gas)
%              PVTG    - properties of wet gas (with vaporised oil)
%              PVDG    - properties of dry gas (no vaporised oil)
%              PVTW    - properties of water
%              DENSITY - fluid density at surface conditions
%              ROCKTAB - rock compaction data
%              ROCK    - rock compressibility
%              SWOF    - gas/oil saturation functions vs gas saturation
%              SGOF    - water/oil saturation functions vs water saturation
%
%              These fields may be added in the future.
%              SWFN    - water saturation functions
%              SOF2    - oil saturation functions (2-phase)
%              SOF3    - oil saturation functions (3-phase)
%              SGFN    - gas saturation functions
%
% COMMENTS:
%
% SEE ALSO:
%   pvt.

%{
#COPYRIGHT#
%}

% $Id: readpvt.m 2338 2009-06-05 17:19:30Z bska $

  opt = struct('verbose', false);
  opt = merge_options(opt, varargin{:});

  action = struct('PVDO',    @(str) readTable(str,3), ...
                  'PVTO',    @readPVTspecial,         ...
                  'PVTW',    @(str) readTable(str,4), ...
                  'PVDG',    @(str) readTable(str,3), ...
                  'PVTG',    @readPVTspecial,         ...
                  'DENSITY', @(str) readTable(str,3), ...
                  'ROCKTAB', @(str) readTable(str,3), ...
                  'ROCK',    @(str) readTable(str,2), ...
                  'SWOF',    @(str) readTable(str,4), ...
                  'SGOF',    @(str) readTable(str,4));%, ...
%                  'SWFN',    @(str) readTable(str, 3),...
%                 'SOF2',    @(str) readTable(str, 3));


  % Scan whole file into string, ignoring comment lines starting with '--';
  [fp, msg] = fopen(filename);
  if fp == -1, error(id('Fopen'), msg), end

  if false

    % Smarter way: not used
    A = textscan(fp,'%s','delimiter','\n','whitespace','');
    m = cellfun(@(x)~isempty(x),regexp(A{1},'^--'));
    S=A{1}(~m);

  else
    S=[];
    while ~feof(fp),  Line = fgets(fp);
      if ~any(regexp(Line, '^--')), S=[S,Line]; %#ok
      end
    end

  end
  fclose(fp);

  % include whitespace in regexp to delimit keywords properly
  keys     = regexp (S, '[A-Z]{2,9}\s','match');
  pos      =  cell2mat(regexp(S, {keys{:}}))';

  % strip whitespace
  keys     = regexprep (keys, '\s+','')';

  [pos, I] = sort(pos, 1); %#ok
  keys     = keys(I);
  pos      = [pos; length(S)+1];

  % For each section identified with keyword, read table
  % into PVT struct
  for i=1:size(keys,1),
    start = pos(i)+length(keys{i});
    stop  = pos(i+1)-1;

    dispif(opt.verbose, 'In %s: reading section %s\n', filename, keys{i});

    try
      Table.(lower(keys{i}))=action.(keys{i})(S(start:stop));
    catch %#ok
      fprintf(2, 'Warning: Keyword %s not read\n', keys{i});
    end
  end

end


function Vector = readNumbers (string)
%-----------------------------------------------------------------------
% Read numbers between '/'

Vector = sscanf(string, '%f');
end


function table = readTable (String, columns)
%-----------------------------------------------------------------------

% Sections
p = [0, regexp(String, '/'), length(String)];


table = [];
for j=1:length(p)-1,
  Num=readNumbers(String(p(j)+1:p(j+1)-1));
  if ~isempty(Num)
    rows      = length(Num)/columns;
    table{j}  = reshape(Num, columns, rows)'; %#ok
    if rows == 1
      % Interpolation routine needs two data points.
      % We interpret a single row as constant data
      table{j}  = [table{j};table{j}]; %#ok
      table{j}(2,1) = table{j} (1,1) + 1;%#ok
    end
  end
end
if numel(table) == 1, table = table{1};end
end



function PVT = readPVTspecial (String)
%-----------------------------------------------------------------------
% Specialized table format used with PVTO keyword.
%
% Read PVTO section with four columns, and two '/' markers between
% Rs-sections and complete tables.
%

%global dostop;
% Strip keyword from string
String = String(length('PVTO')+1:end);

% Sections
Pstart = [regexp(String, '/\s*/','start'), length(String)+1];
Pend   = [1,regexp(String, '/\s*/','end')];


for i=1:length(Pstart)-1,

  S = String(Pend(i):Pstart(i)-1);

  % Subsections
  p = [0, regexp(S, '/'), length(S)+1];

  Tab = [];
  for j=1:length(p)-1,

    Num=readNumbers(S(p(j)+1:p(j+1)-1));

    if ~isempty(Num)
      rows = (length(Num)-1)/3;
      R =Num(1);
      Tab=[Tab; [[R ;NaN(rows-1, 1)],reshape(Num(2:end), 3, rows)']];%#ok
    end
  end
  PVT{i} = Tab;%#ok
end
if numel(PVT) == 1, PVT=PVT{1}; end
end

%--------------------------------------------------------------------------

function s = id(s)
s = ['readpvt:', s];
end
