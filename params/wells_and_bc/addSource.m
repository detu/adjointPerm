function src = addSource(src, c, r, varargin)
%Add an explicit source to (new or existing) source object.
%
% SYNOPSIS:
%   src = addSource(src, cells, values)
%   src = addSource(src, cells, values, 'pn1', pv1)
%
% PARAMETERS:
%   src    - Source structure from a prior call to 'addSource' which will
%            be updated on output or an empty array (src==[]) in which case
%            a new source structure is created.
%
%   cells  - Indices in external model for which this source should be
%            applied.
%
%   values - Strength of source.  One scalar value for each cell in
%            'cells'.  Note that the values in 'values' are interpreted as
%            flux rates (typically in units of m^3/Day) rather than as flux
%            density rates (which must be integrated over the cell volumes
%            in order to obtain flux rates).  Specifically, the mimetic
%            pressure solvers do not integrate these values.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   sat    - Fluid composition of injected fluid in injection cells.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of cells in 'cells' and 'm' the number of components in
%            the saturation. For m=3, the columns interpreted as
%              1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            'blackoilUpwFE' and will be ignored for production source
%            cells (i.e. when values < 0).
%
%            Default value: sat = [] (assume single-phase flow).
%
% RETURNS:
%   src - New or updated source structure having the following fields:
%           - cell -- Cells for which explicit sources are provided
%           - rate -- Rates or values of these explicit sources
%           - sat  -- Fluid composition of injected fluids in injection
%                     cells.
%
% EXAMPLE:
%   See simpleSRCandBC.m
%
% SEE ALSO:
%   addWell, addBC, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Id: addSource.m 2338 2009-06-05 17:19:30Z bska $

opt = struct('sat', []);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

if isempty(src),
   src = struct('cell', [], 'rate', [], 'sat', []);
end

assert (numel(c) == numel(r));
assert (numel(s)==0 || (size(s,2) >= 1 && numel(r) == size(s,1)));
assert ((size(s,2) == size(src.sat,2)) || (size(src.sat,2) == 0));

src.cell = [src.cell; c(:)];
src.rate = [src.rate; r(:)];
src.sat  = [src.sat ; s];
