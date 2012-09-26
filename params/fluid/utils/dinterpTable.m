function dyidx = dinterpTable(x, y, xi, varargin)
%Compute derivative of one-dimensional interpolant, possibly using splines.
%
% SYNOPSIS:
%   dyi = dinterpTable(X, Y, xi)
%   dyi = dinterpTable(X, Y, xi, 'pn1', pv1, ...)
%
% PARAMETERS:
%   X       - Nodes at which underlying function y=y(x) is sampled.
%
%   Y       - Values of the underlying function y=y(x) at the nodes, X.
%
%   xi      - Evaluation points for new, interpolated, values of the
%             derivative y'(x).
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%                - spline -- Whether or not to use spline interpolation.
%                            Logical.  Default value: spline=false (use
%                            linear interpolation/extrapolation).
%
% RETURNS:
%   dyi - Approximate (interpolated/extrapolated) values of the derivative
%         of the function y=y(x) at the points xi.
%
% SEE ALSO:
%   dinterpq1, interpTable, interp1, spline.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

   opt = struct('spline', false);
   opt = merge_options(opt, varargin{:});

   if opt.spline
      der       = spline(x, y);
      d         = diag(der.order-1 : -1 : 1, 1);
      der.coefs = der.coefs * d;
      dyidx     = ppval(der, xi);
   else
      dyidx     = dinterpq1(x, y, xi);
   end
end
