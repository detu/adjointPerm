function [krg, kro] = sgof(T, varargin)
%Construct gas/oil relperm evaluation functions from SGOF table.
%
% SYNOPSIS:
%   [krg, kro] = sgof(table)
%   [krg, kro] = sgof(table, Swco)
%
% PARAMETERS:
%   table     - Black Oil PVT/Relperm table structure as defined by
%               function 'readpvt'.  Must contain valid field 'sgof'.
%
%   Swco      - Connate water saturation for runs containing water as an
%               active phase.  Maximum gas saturation in table and
%               simulation run must not exceed Sg_{max} = 1 - Swco.
%               OPTIONAL.  Default value: Swco = 0 (no water).
%
% RETURNS:
%   krg - Two-component cell array of gas relative and differentiated
%         relative permeability evaluation functions as single-variate
%         functions of the gas saturation (Sg).
%
%   kro - Two-component cell array of oil relative and differentiated
%         relative permeability evaluation functions as single-variate
%         functions of the oil saturation (So).
%
% SEE ALSO:
%   readpvt, swof, initBlackOilRelPerm.

% $Date: 2009-06-10 19:08:19 +0200 (on, 10 jun 2009) $
% $Revision: 2357 $

   if ~isfield(T, 'sgof'),
      error(msgid('SGOF:invalid'), ...
            'Relative permeability table must have valid field ''sgof''.');
   end

   if iscell(T.sgof),
      if numel(T.sgof) > 1,
         error(msgid('SGOF:TooMany'), ...
               '''sgof'' field must contain a single SGOF table.');
      end
      T.sgof = T.sgof{1};
   end

   % Set connate water saturation (if water is present).
   Swco = 0;
   if nargin > 1 && isnumeric(varargin{1}) && numel(varargin{1}) == 1,
      Swco = varargin{1};
   end
   max_sat = max(T.sgof(:,1));

   if     max_sat > 1 - Swco,
      warning(msgid('MaxSat:NoSwco'), ...
             ['Maximum gas saturation in ''SGOF'' table (=%f) does ', ...
              'not account for presence of connate water.'], max_sat);
   elseif max_sat < 1 - Swco,
      warning(msgid('MaxSat:LessSwco'), ...
             ['Maximum gas saturation in ''SGOF'' table (=%f) is ',   ...
              'less than connate water saturation.\n',                ...
              'The maximum gas saturation should normally be 1-Swco', ...
              ' (=%f).'], max_sat, 1 - Swco);
   end

   % Assert Swco \in [0, 1].
   assert (~((Swco < 0) || (Swco > 1)));

   % Basic validation of input data.
   assert (~(abs(T.sgof(1,1)) > 0));      % min(Sg) == 0.
   assert (all(diff(T.sgof(:,1)) > 0));   % Sg monotonically increasing.

   assert (~(abs(T.sgof( 1 ,2)) > 0));    % krg(min(Sg))  == 0.
   assert (~(abs(T.sgof(end,3)) > 0));    % kro(Sg_{max}) == 0.

   assert (all (diff(T.sgof(:,1)) > 0));  % Sg monotonically increasing.
   assert (~any(diff(T.sgof(:,2)) < 0));  % Level or increasing down column
   assert (~any(diff(T.sgof(:,3)) > 0));  % Level or decreasing down column

   % Gas relative permeability as a function of water saturation, sg.
   % Note: Gas saturation is always <= 1-Swco.
   function varargout = interp_gas(sg)
      kr = {(@(sg)  interpTable(T.sgof(:,1), T.sgof(:,2), ...
                                min(sg, 1 - Swco))), ...
            (@(sg) dinterpTable(T.sgof(:,1), T.sgof(:,2), ...
                                min(sg, 1 - Swco)))};

      varargout = cellfun(@(fun) fun(sg), { kr{1 : nargout} }, ...
                          'UniformOutput', false);
   end
   krg = @interp_gas;

   % Oil relative permeability as a function of oil saturation, so.
   % Note: Oil saturation is always <= 1-Swco.
   function varargout = interp_oil(so)
      kr = {(@(so)   interpTable(T.sgof(:,1), T.sgof(:,3), ...
                                 max(1 - Swco - so, 0))),       ...
            (@(so) -dinterpTable(T.sgof(:,1), T.sgof(:,3), ...
                                 max(1 - Swco - so, 0)))};

      varargout = cellfun(@(fun) fun(so), { kr{1 : nargout} }, ...
                          'UniformOutput', false);
   end
   kro = @interp_oil;
end
