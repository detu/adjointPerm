function t = dynamicCoarseWeight(cellno, BI, ct, phi, v, dpdt)
%Compute synthetic multiscale weighting function.
%
% SYNOPSIS:
%   theta = dynamicCoarseWeight(cellno, BI, ct, phi, v, dpdt)
%
% PARAMETERS:
%   cellno - A map giving the global cell number of any given half-face.
%            That is, cellno(i) is the global cell to which half-face 'i'
%            is connected.  Typically,
%
%                nc     = G.cells.num;
%                nf     = double(G.cells.numFaces);
%                cellno = rldecode((1 : nc) .', nf);
%
%            when 'G' is the grid of a reservoir model.
%
%   BI     - Inverse mass matrix, correctly updated for effects of total
%            mobility, such that the expression
%
%                f' * (BI \ f)
%
%            gives the (squared) energy norm of a function 'f' represented
%            on all half faces of the model.
%
%   ct     - Total compressibility.  One scalar value for each global cell
%            in the reservoir model.
%
%   phi    - Porosity.  One scalar value for each global cell in the model.
%
%   v      - Current half-face fluxes for all cells in the model.
%
%   dpdt   - An approximate time derivative of the cell pressure values.
%            One scalar value for each cell in the grid.
%
% RETURNS:
%   theta - A \theta function suitable for passing in option pair
%           ('BasisWeighting',theta) to function generateCoarseSystem.
%
% SEE ALSO:
%   computeMimeticIP, generateCoarseSystem, pvt, rldecode.

%{
#COPYRIGHT#
%}

% $Id: dynamicCoarseWeight.m 1953 2009-03-31 10:54:12Z bska $

t = phi .* ct .* dpdt;
t = accumarray(cellno, v .* (BI \ v)) - t;
