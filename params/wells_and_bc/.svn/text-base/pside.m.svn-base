function bc = pside(bc, G, side, I1, I2, pressure, varargin)
%Impose pressure boundary condition on global side.
%
% SYNOPSIS:
%   bc = pside(bc, G, side, I1, I2, p)
%   bc = pside(bc, G, side, I1, I2, p, 'pn', pv)
%
% PARAMETERS:
%   bc     - boundary condition structure as defined by function 'addBC'.
%
%   G      - Grid structure as described by grid_structure.  Currently
%            restricted to grids produced by functions cartGrid and
%            tensorGrid.
%
%   side   - Global side from which to extract face indices.
%            Must be one of {'LEFT', 'RIGHT', 'FRONT', 'BACK', ...
%                            'BOTTOM', 'TOP'}
%
%   I1, I2 - Cell index ranges for local axes one and two, respectively.
%            Must contain *ALL* indices (e.g., I1 = 1:50, I2 = 1:4).
%
%   p      - Pressure value, in units of Pascal, to be applied to the face.
%            Either a scalar or a vector of NUMEL(I1)*NUMEL(I2) values.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   sat    - Fluid composition of fluid injected across inflow faces.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of individual faces specified by (I1,I2) (i.e.,
%            n==NUMEL(I1)*NUMEL(I2)) or one.  If m=3, the columns of 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            'blackoilUpwFE' and will be ignored for outflow faces.
%
%            Default value: sat = [] (assume single-phase flow).
%
%   range  - Search range in the perpendicular direction. Restricts the
%            search of outer faces to a subset of the cells in the model.
%            Default value: range = [] (do not restrict search)
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLE:
%   See simpleBC, simpleSRCandBC.
%
% SEE ALSO:
%   fluxside, addBC, solveIncompFlow, grid_structure.

%{
#COPYRIGHT#
%}

% $Id: pside.m 2338 2009-06-05 17:19:30Z bska $

error(nargchk(6, 10, nargin, 'struct'))

opt = struct('sat', [], 'range', []);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = boundaryFaceIndices(G, side, I1, I2, opt.range);

assert (any(numel(pressure) == [1, numel(ix)]));
assert (numel(sat)==0 || any(size(sat,1) == [1, numel(ix)]));

if size(sat,1)     == 1, sat      = sat(ones([numel(ix), 1]), :); end
if numel(pressure) == 1, pressure = pressure(ones([numel(ix), 1])); end

bc = addBC(bc, ix, 'pressure', pressure, 'sat', sat);
