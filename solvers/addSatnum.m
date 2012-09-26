function rSol = addSatnum(rSol, G, grdecl)
% Add satnum fields to reservoir solution data structure.
%
% SYNOPSIS:
%   resSol = addSatnum(resSol, G, grdecl)
%
% DESCRIPTION:
%   Adds fields satnum and satnumActive to the reservoir solution
%   structure.
%
% PARAMETERS:
%   resSol  - Reservoir soultion structure created by initResSol 
%   G       - Grid data structure.
%   grdecl  - Raw pillar grid structure, as defined by function
%            'readGRDECL', must have fields ACTNUM and SATNUM.
% 
% RETURNS:
%   resSol - Updated reservoir solution structure with new fields: 
%              - satnum       -- logical matrix of size
%                                G.cells.num-by-satnumActive
%
%              - satnumActive -- Array of 
%                                active satnum = unique(grdecl.SATNUM) 
%
% SEE ALSO:
%   initWellSol, initResSol, initSatnumFluid

%{
#COPYRIGHT#
%}

% $Id: addSatnum.m $

assert ( all(isfield(grdecl, {'SATNUM', 'ACTNUM'})))

satnum = grdecl.SATNUM(G.cells.indexMap);
satnum = satnum(logical(grdecl.ACTNUM(G.cells.indexMap)));

% identify unique satnum values
active = unique(satnum);

inx         = zeros(max(satnum),1);
inx(active) = 1:numel(active);

% make index matrix of satnum values in grid
rSol.satnum = logical(sparse(1:G.cells.num, inx(satnum), 1, ...
                      G.cells.num, numel(active)));

rSol.satnumActive = active; 
