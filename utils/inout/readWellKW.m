function w  = readWellKW(fid, w, kw)
%Read well definitions from an ECLIPSE Deck specification.
%
% SYNOPSIS:
%   w = readWellKW(fid, w, kw)
%
% PARAMETERS:
%   fid - File identifier
%   w   - Struct array with well definition.
%   kw  - Keyword
%
% RETURNS:
%   w   - Updated struct array with well definition.
%
% COMMENTS:
%         The currently recognized keywords are:
%         'WELSPECS', 'COMPDAT', 'WCONINJE',  'WCONPROD'
%
% SEE ALSO:
%   readGRDECL, processWells.

%{
#COPYRIGHT#
%}

% $Id: readWellKW.m 2646 2009-09-04 09:32:08Z hnil $

%% Define the records to be read for each kw.
%
switch kw                       
   case 'WELSPECS'
      data = {'I', 'J'};
      number = [3 4];
      numeric= [1 1];
   case 'COMPDAT'
      data = { 'I_vec', 'J_vec', 'K_vec_u', 'K_vec_l', ...
               'state_vec', 'transm_vec', 'diam_vec', 'Kh_vec', ...
               'skin_vec', 'dir_vec'};
      number = [2 3 4 5 6 8 9 10 11 13];
      numeric= [1 1 1 1 0 1 1 1 1 0];
   case 'WCONINJE'
      data = { 'inj_type', 'state', 'control', 'surf_vol_rate', ...
               'vol_rate', 'bhp'};
      number = [2 3 4 5 6 7];
      numeric= [0 0 0 1 1 1];
      well_type = 'Injection';
   case 'WCONPROD'
      data = { 'state', 'control', 'oil_rate', 'water_rate', 'gas_rate', ...
               'liquid_rate', 'vol_rate', 'bhp'};
      number = [2 3 4 5 6 7 8 9];
      numeric= [0 0 1 1 1 1 1 1];
      well_type = 'Production';
    case 'WCONHIST'
      data = { 'state', 'control', 'oil_rate', 'water_rate', 'gas_rate'};
      number = [2 3 4 5 6 7 8 9];
      numeric= [0 0 1 1 1 1 1 1];
      well_type = 'Production';
      
end

% Read first line.
lin = fgetl(fid);

% Main loop. -------------------------------------------------------------- 
% Loop lines until / (end of record).
%
while isempty(regexp(lin, '^/', 'once'))
   count = 1;
   % Skip blank lines and comments.
   if isempty(lin) || ~isempty(regexp(lin, '^--', 'match'))
      lin = fgetl(fid);
      continue
   end

   % Split line and remove ' '
   split = regexp(lin, '(\w\.*\-*\**)+', 'match');
   
   
   %% Process first part of 'lin'. ---------------------------------------- 
   %
   wellName = split(1);

   % Check if well is defined by searching well names in struct.
   % k = index to well in array of structs or empty if new well.
   names = {w.name};
   k = strmatch(wellName, names);
   
   % Initialize well from default values if well is not already defined.
   if isempty(k)
      k= size(w, 2) + 1;
      w(k) = w(1);
      w(k).name = char(wellName);
      if ~strcmp(kw, 'WELSPECS')
         error('readWellKW:Input:Unknown', ...
              ['Well `%s'' referenced from keyword %s before it is', ...
               'defined with keyword WELSPECS'], w(k).name, kw);
      end
      
   % Add defaults for fields in COMPDAT.
   elseif strcmp(kw, 'COMPDAT') 
     if w(k).compdat_count > 0, 
         w(k) = addDefaults(w(k), w(1), data, kw);
         count = w(k).compdat_count;
     else w(k).compdat_count = 1; % first time defaults are ok.   
     end        
   
   % Define well type and check if prod/inj data already is added for well.   
   elseif (strcmp(kw, 'WCONINJE') || strcmp(kw, 'WCONPROD') ||  strcmp(kw, 'WCONHIST'))
     if ~isempty(w(k).type)
       warning('readWellKW:Input:Unknown', ...
              ['Production or injection data from Well `%s'' is', ...
               ' already defined \n' 'but will be deleted and replaced',...
               'by values: \n %s'], w(k).name, lin);
            w(k) = addDefaults(w(k), w(1), data, kw);
       
     end      
      w(k).type = well_type;
   end
   
   % Define indices for looping remaing entries of 'lin'
   i = 2;   % Index to data in array split.
   j = 2;   % Index to data in a full keyword record.
   if strcmp(kw, 'WELSPECS')  && ...
         isempty(regexp(split(2), '[A-Za-z]', 'once'))
      % Inidices must be updated when field name is empty.
      i = i + 1;
      j = j + 1;
   end

   %% Process remaining entries of 'lin' ----------------------------------
   %
   while i < length(split) % numbers
      next = char(split(i));
      ind  = regexp(next, '[*]' );
      
      % Next is on the form n*, defaulting the next n values.
      if ind 
         n = str2num(next(1:ind-1)); %#ok
         i = i + 1;
         j = j + n;
         % Continue to next entry in split.   
         
      else % Next is a value
         % inx = index to data in the arrays: data, number and numeric.
         inx = find(number == j);
         if ~isempty(inx) && (i < length(split))
            % Next is a value to be put in the well-structure.
            val = next; 
            if numeric(inx), 
               % Convert char to num.
               val = str2num(val);  %#ok 
            end 
            
            % Read a value into a vector.
            if strfind(data{inx}, 'vec')
               % Convert 'OPEN'/'CLOSED' to numerical val.
               if strcmpi(val, 'OPEN'), val = 1;
               elseif strcmpi(val, 'CLOSED'), val = 0; 
               end
               w(k).(data{inx})(count) = val;
            else
               w(k).(data{inx}) = val;
            end
         end
         i = i + 1;
         j = j + 1;
      end
   end
   % Calculate values where needed
   if(strcmp(kw, 'WCONHIST'))
         w(k).liquid_rate = w(k).oil_rate +  w(k).water_rate;
         w(k).vol_rate = w(k).liquid_rate + w(k).gas_rate;
   end
   % Read next line.
   lin = fgetl(fid);
end
end

function w = addDefaults(w, def_w, data, kw)
% Put defaults in vector fields that belong to keyword kw.

if strcmp(kw, 'COMPDAT')
   % Put defaults and update count value. 
   i = w.compdat_count + 1;
   w.compdat_count = i;
   for k = 1: numel(data)
      w.(data{k})(i) = def_w.(data{k});
   end

else % Put defaults.
   for k = 1: numel(data)
      w.(data{k}) = def_w.(data{k});
   end
end
end





