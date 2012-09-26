function resSol = initResSol(G, p0, s0, z0)
%Initialize reservoir solution data structure.
%
% SYNOPSIS:
%   resSol = initResSol(G, p0, s0, z0)
%
% DESCRIPTION:
%   Initialize the reservoir solution structure to uniform cell and face
%   pressures and no-flow flux fields throughout the reservoir.  Also, set
%   all-zero phase saturations and phase densities in all cells in model.
%
%
% PARAMETERS:
%   G  - Grid data structure.
%   p0 - Initial uniform reservoir pressure (scalar).
%   s0 - Initial uniform reservoir saturation (default: s0=0)
%
% RETURNS:
%   resSol - Initialized reservoir solution structure having fields
%              - cellPressure -- REPMAT(p0, [G.cells.num, 1])
%              - facePressure -- REPMAT(p0, [G.faces.num, 1])
%              - cellFlux     -- Initial, all-zero cell fluxes for all
%                                cells in model.
%              - faceFlux     -- Initial, all-zero face fluxes for all
%                                faces in model.
%              - s            -- REPMAT(s0, [G.cells.num, 1])
%              - z            -- REPMAT(s0, [G.cells.num, 1])
%
% REMARKS:
%   resSol.s   - An n-by-m array of fluid compositions with 'n' being the
%                number of faces in 'faces' and for m=3, the columns
%                interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor. This
%                field is for the benefit of transport solvers.
%
%   The field resSol.z is only added if the saturation has three or more
%   components. 
%
% SEE ALSO:
%   initWellSol, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Id: initResSol.m 2338 2009-06-05 17:19:30Z bska $

if nargin <= 2,    s0 = 0; 
end
   
[nc, nf] = deal(G.cells.num, G.faces.num);

s0 = reshape(s0, 1, []);

resSol = struct('cellPressure', repmat(p0, [nc, 1]),              ...
                'facePressure', repmat(p0, [nf, 1]),              ...
                'cellFlux',     zeros([size(G.cellFaces,1), 1]),  ...
                'faceFlux',     zeros([nf, 1]),                   ...
                's',            s0(ones([nc, 1]), :));
if nargin == 4, resSol.z = z0(ones([nc, 1]), :); end
