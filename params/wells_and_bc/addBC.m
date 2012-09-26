function bc = addBC(bc, f, t, v, varargin)
%Add boundary condition to (new or existing) BC object
%
% SYNOPSIS:
%   bc = addBC(bc, faces, type, values)
%   bc = addBC(bc, faces, type, values, 'pn1', pv1, ...)
%
% PARAMETERS:
%   bc     - Boundary condition structure from a prior call to 'addBC'
%            which will be updated on output or an empty array (bc==[]) in
%            which case a new boundary condition structure is created.
%
%   faces  - Global faces in external model for which this boundary
%            condition should be applied.
%
%   type   - Type of boundary condition.  Supported values are 'pressure'
%            and 'flux'.
%
%   values - Boundary condition value.  Interpreted as a pressure value (in
%            units of 'Pa') when type=='pressure' and as a flux value (in
%            units of 'm^3/s') when type=='flux'.  One scalar value for
%            each face in 'faces'.
%
%            Note: If type=='flux', the values in 'values' are interpreted
%            as injection fluxes (into the reservoir).  Specifically, a
%            positive value is interpreted as an injection flux.  To
%            specify an extraction flux (i.e., flux out of the reservoir),
%            the caller should provide a negative value in 'values'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   sat    - Fluid composition of fluid injected across inflow faces.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of faces in 'faces' and for m=3, the columns
%            interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            'blackoilUpwFE' and will be ignored for outflow faces.
%
%            Default value: sat = [] (assume single-phase flow).
%
% NOTE
%   For convenience, values and sat may contain a single value.  This value
%   is then used for all faces specified in the call.
%
% RETURNS:
%   bc - New or updated boundary condition structure having the following
%        fields:
%           - face  -- External faces for which explicit BCs are provided.
%           - type  -- Cell array of strings denoting type of BC.
%           - value -- Boundary condition values for all faces in 'face'.
%           - sat   -- Fluid composition of fluids passing through inflow
%                      faces.
%
% SEE ALSO:
%   pside, fluxside, addSource, addWell, solveIncompFlow, grid_structure.

%{
#COPYRIGHT#
%}

% $Id: addBC.m 2975 2009-10-12 08:08:04Z ilig $

opt = struct('sat', []);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

if isempty(bc),
   bc = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
end

% check that type is valid
assert (strcmp(t, 'pressure') || strcmp(t, 'flux'));
% check that saturation is valid
assert (numel(s)==0 || size(s,2) >= 1);
assert ((size(s,2) == size(bc.sat,2)) || (size(bc.sat,2) == 0));

t = { t };

if size(s,1) == 1, s = s(ones([numel(f), 1]), :); end
if numel(v)  == 1, v = v(ones([numel(f), 1]));    end

% check that v and s are same length as faces (or s empty)
assert(numel(v) == numel(f));
assert(any(size(s,1) == [0, numel(f)]));

bc.face  = [bc.face   ; f(:)];
bc.type  = {bc.type{:}, t{ones([1, numel(f)])}};
bc.value = [bc.value  ; v(:)];
bc.sat   = [bc.sat    ; s];
