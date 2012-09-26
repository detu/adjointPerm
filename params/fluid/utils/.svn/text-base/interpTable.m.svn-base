function yi = interpTable(X, Y, xi, varargin)
%Interpolate a one-dimensional table, possibly using splines.
%
% SYNOPSIS:
%   yi = interpTable(X, Y, xi)
%   yi = interpTable(X, Y, xi, 'pn1', pv1, ...)
%
% PARAMETERS:
%   X       - Nodes at which underlying function y=y(x) is sampled.
%
%   Y       - Values of the underlying function y=y(x) at the nodes, X.
%
%   xi      - Evaluation points for new, interpolated, function values.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%                - spline -- Whether or not to use spline interpolation.
%                            Logical.  Default value: spline=false (use
%                            linear interpolation/extrapolation).
%
% RETURNS:
%   yi - Approximate (interpolated/extrapolated) values of the function
%        y=y(x) at the points xi.
%
% SEE ALSO:
%   dinterpTable, interp1, spline.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

   opt = struct('spline', false);
   opt = merge_options(opt, varargin{:});

   if opt.spline,
      method = 'spline';
   else
      method = 'linear';
   end

   yi = interp1(X, Y, xi, method, 'extrap');
end
