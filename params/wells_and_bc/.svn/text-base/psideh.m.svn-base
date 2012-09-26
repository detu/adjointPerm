function bc = psideh(bc, G, side, I1, I2, fluid, varargin)
%Impose Dirichlet boundary condition (pressure) on global side.
%
% SYNOPSIS:
%   bc = psideh(bc, G, side, I1, I2, p)
%   bc = psideh(bc, G, side, I1, I2, p, 'pn', pv)
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
%   sat    - Fluid composition of fluid outside of the reservoir.
%            An m array of fluid phase saturations. If m=3, 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is to side the density of the outside fluid and to
%            set the saturation of incoming fluid in a transport solver
%
%            Default value: sat = 0 (assume single-phase flow).
%
%   range  - Search range in the perpendicular direction. Restricts the
%            search of outer faces to a subset of the cells in the model.
%            Default value: range = [] (do not restrict search)
%
%   ref_depth - Refrance depth for pressure. Default i 0
%
%   ref_pressure - Referance pressure. Default is 0
%
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLE:
%   See simpleBC.m, simpleSRCandBC.m.
%
% SEE ALSO:
%   pside, fluxside, addBC, solveIncompFlow, grid_structure.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

error(nargchk(6, 10, nargin, 'struct'))

opt = struct('sat', 0, 'range', [],'ref_position',[0.0,0.0,0.0],'ref_pressure',0.0);
opt = merge_options(opt, varargin{:});
sat = opt.sat;

ix = boundaryFaceIndices(G, side, I1, I2, opt.range);

%assert (any(numel(pressure) == [1, numel(ix)]));
assert (numel(sat)==0 || any(size(sat,1) == [1, numel(ix)]));
assert( size(sat,1) == 1)
if size(sat,1)     == 1, sat      = sat(ones([numel(ix), 1]), :); end
%if numel(pressure) == 1, pressure = pressure(ones([numel(ix), 1])); end
%keyboard
pressure = sum( (G.faces.centroids(ix,:)-repmat(opt.ref_position,length(ix),1)).*repmat(gravity(),length(ix),1) ,2).*fluid.omega(struct('s',sat))+opt.ref_pressure;

bc = addBC(bc, ix, 'pressure', pressure, 'sat', sat);
