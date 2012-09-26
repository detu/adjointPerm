function pv = poreVolume(G, rock, varargin)
%Compute pore volumes of individual cells in grid.
%
% SYNOPSIS:
%   pv = poreVolume(G, rock)
%
% PARAMETER:
%   G    - Grid data structure.
%          Must contain valid field 'G.cells.volumes'.
%
%   rock - Rock data structure.
%          Must contain valid field 'rock.poro'.
%
% RETURNS:
%   pv   - Vector of size (G.cells.num)-by-1 of pore volumes for each
%          individual cell in the grid.  This typically amounts to the
%          expression rock.poro .* G.cells.volumes, but the function
%          handles non-unit net-to-gross factors as well.
%
% SEE ALSO:
%  computeGeometry.

%{
#COPYRIGHT#
%}

% $Id: poreVolume.m 1951 2009-03-31 10:28:32Z bska $

pv = rock.poro .* G.cells.volumes;

if isfield(rock, 'ntg'),
   pv = pv .* rock.ntg;
end
