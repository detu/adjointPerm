function fluid = initBlackoilFluid(T, varargin)
%Return function that evaluate Black Oil PVT model
%
% SYNOPSIS:
%   fun = initBlackoilFluid(T, varargin)
%
% PARAMETERS:
%   T       - Struct of (Eclipse) PVT tables.  It must contain
%             .pvto(.pvdo), .pvtg(.pvdg), .pvtw, .density fields.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - verbose -- Whether or not to emit informational messages
%                            while interpolating PVT tables.
%                            Logical.  Default value = FALSE.
%
% RETURNS:
%   fluid.pvt- Function that evaluates pvt data for the Black Oil.
%             Given n-by-1 pressures p and n-by-np surface volumes z,
%             the call
%
%                  [c, rho, mu, S] = fluid.pvt(p, z)
%
%             computes respectively n-by-np phase compressibilities (c),
%             densities (rho), viscosities (mu) and saturations (S).
%
%   fluid.relperm - Function that evaluates relative permeability curves
%             specified by SWOF and SGOF.
%
% COMMENTS:
%   The individual Black Oil phases in these arrays are ordered into
%   columns as follows:
%
%     'Aqua'   -> Column 1
%     'Liquid' -> Column 2
%     'Vapor'  -> Column 3
%
%   For example, B(i,2) is the volume formation factor of the 'Liquid'
%   phase in cell 'i' while c(i,3) is the 'Vapor' phase compressibility in
%   cell 'i'.
%
% EXAMPLE:
%
% SEE ALSO:
%   initBlackOilPVT, initBlackOilRelPerm.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
% $Revision: 2926 $

   [pvt, rhoS, miscible, info] = initBlackOilPVT    (T, varargin{:});
   fluid.relperm               = initBlackOilRelPerm(T, varargin{:});

   fluid.pvt                   = pvt;
   fluid.surfaceDensity        = rhoS;
   fluid.miscible              = miscible;
   fluid.info                  = info;
end
