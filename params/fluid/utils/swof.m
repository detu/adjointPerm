function [krw, kro, Swco, krocw] = swof(T, varargin)
%Construct water/oil relperm evaluation functions from SWOF table.
%
% SYNOPSIS:
%   [krw, kro, Swco, krocw] = swof(table)
%
% PARAMETERS:
%   table - Black Oil PVT/Relperm table structure as defined by function
%           'readpvt'.  Must contain valid field 'swof'.
%
% RETURNS:
%   krw   - Two-component cell array of Water relative and differentiated
%           relative permeability evaluation functions as single-variate
%           functions of the water saturation (Sw).
%
%   kro   - Two-component cell array of oil relative and differentiated
%           relative permeability evaluation functions as single-variate
%           functions of the oil saturation (So).
%
%   Swco  - Connate water saturation.  Oil saturation in run must not
%           exceed 1-Swco.
%
%   krocw - Oil relative permeability at So=1-Swco (at connate water).
%
% NOTE:
%   Function 'swof' does presently only support a single saturation
%   function througout the reservoir.  In other words the SATNUM keyword is
%   not presently recognised.
%
% SEE ALSO:
%   readpvt, sgof, initBlackOilRelPerm.

% $Date: 2009-08-31 17:14:49 +0200 (ma, 31 aug 2009) $
% $Revision: 2613 $

   if ~isfield(T, 'swof'),
      error(msgid('SWOF:invalid'), ...
            'Relative permeability table must have valid field ''swof''.');
   end

   if iscell(T.swof),
      if numel(T.swof) > 1,
         error(msgid('SWOF:TooMany'), ...
               '''swof'' field must contain a single SWOF table.');
      end
      T.swof = T.swof{1};
   end

   % Basic validation of input data.
   assert (~(abs(T.swof( 1 ,2)) > 0));    % krw(Swco)     == 0
   assert (~(abs(T.swof(end,3)) > 0));    % kro(Sw_{max}) == 0

   assert (all (diff(T.swof(:,1)) > 0));  % Sw monotonically increasing
   assert (~any(diff(T.swof(:,2)) < 0));  % Level or increasing down column
   assert (~any(diff(T.swof(:,3)) > 0));  % Level or decreasing down column

   % Connate water saturation (>= 0) is first Sw encountered in table.
   Swco = T.swof(1,1);  assert (~(Swco < 0));

   % Oil relperm at Sw = Swco (=> So = 1 - Swco)
   krocw = T.swof(1,3);

   if Swco > 0, T.swof = [0, T.swof(1,2:end); T.swof]; end
   
   % Water relative permeability as a function of water saturation, sw.
   % Note: Water saturation is always >= Swco.
   function varargout = interp_water(sw)
      kr = {(@(sw)  interpTable(T.swof(:,1), T.swof(:,2), ...
                                max(sw, Swco))), ...
            (@(sw) dinterpTable(T.swof(:,1), T.swof(:,2), ...
                                max(sw, Swco)))};

      varargout = cellfun(@(fun) fun(sw), { kr{1:nargout} }, ...
                          'UniformOutput', false);
   end
   krw = @interp_water;

   % Oil relative permeability as a function of oil saturation, so.
   % Note: Oil saturation is always <= 1 - Swco.
   function varargout = interp_oil(so)
      kr = {(@(so)   interpTable(T.swof(:,1), T.swof(:,3), ...
                                 max(1 - Swco - so, 0))),       ...
            (@(so) -dinterpTable(T.swof(:,1), T.swof(:,3), ...
                                 max(1 - Swco - so, 0)))};

      varargout = cellfun(@(fun) fun(so), { kr{1:nargout} }, ...
                          'UniformOutput', false);
   end
   kro = @interp_oil;
end
