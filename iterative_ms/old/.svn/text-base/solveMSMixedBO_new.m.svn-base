function [resSol, wellSol] = solveMSMixedBO(resSol, wellSol, G, CG, p, ...
                                            rock, S, CS, fluid, p0, dt, ...
                                            varargin)
%One iteration of succsubst algorithm on Black Oil pressure system.
%
% SYNOPSIS:
%   [xr, xw] = solveBlackOilWellSystem(xr, xw, G, rock, S, fluid, p0, dt)
%   [xr, xw] = solveBlackOilWellSystem(xr, xw, G, rock, S, fluid, p0, dt, ...
%                                      'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   xr, xw  - Reservoir and well solution structures either properly
%             initialized or the results from a previous call to function
%             'solveBlackOilWellSystem' (i.e., previous iteration of
%             successive substitution algorithm or previous time step),
%             and, possibly, a transport solver.
%
%   G       - Grid structure as described in grid_structure.
%
%   rock    - Rock data structure.  Must contain valid field 'rock.perm'.
%
%   S       - Linear system structure as defined by function
%             'computeMimeticIP'.
%
%   fluid   - Black Oil fluid object as defined by, e.g., function
%             'initBlackoilFluid'.
%
%   p0      - Vector, length G.cells.num, of cell pressures at previous
%             time step (not previous iteration of successive substitution
%             algorithm).
%
%   dt      - Time step size.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells    - Well structure as defined by function 'addWell'.  May be
%              empty (i.e., W = [], default value) which is interpreted as
%              a model without any wells.
%
%   bc       - Boundary condition structure as defined by function 'addBC'.
%              This structure accounts for all external boundary conditions
%              to the reservoir flow.  May be empty (i.e., bc = [], default
%              value) which is interpreted as all external no-flow
%              (homogeneous Neumann) conditions.
%
%   src      - Explicit source contributions as defined by function
%              'addSource'.  May be empty (i.e., src = [], default value)
%              which is interpreted as a reservoir model without explicit
%              sources.
%
%   LinSolve - Handle to linear system solver software to which the fully
%              assembled system of linear equations will be passed.
%              Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%              in order to solve a system Ax=b of linear equations.
%              Default value: LinSolve = @mldivide (backslash).
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - cellPressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%          - cellFlux     -- Outgoing flux across each local interface for
%                            all cells in the model.
%          - faceFlux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% SEE ALSO:
%   initBlackoilFluid, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-07 19:22:40 +0200 (sø, 07 jun 2009) $
% $Revision: 2342 $

opt = struct('bc', [], 'src', [], 'wells', [], ...
    'LinSolve', @mldivide);
opt = merge_options(opt, varargin{:});

W = opt.wells;

[B, C, D, P, f, g, h] = ...
    hybLinSys(resSol, G, S, rock, fluid, p0, dt, ...
              'bc', opt.bc, 'src', opt.src);
[Bw, Cw, Dw, fw, gw, hw] = hybLinWellSys(resSol, wellSol, G, S, W, rock, fluid);

B = blkdiag(B, Bw);
C = [C ; Cw];
[D, Do, h] = mixedMappings(G, W, D, Dw, h, hw, opt.bc);

f = [f; fw]; 
g = g + gw;

[flux, p, lamN] = solveMixedLinSys(B, C, D, P, f, g, h, Do);

[resSol, wellSol] = packSol(G, W, resSol,wellSol, flux, p, lamN, opt.bc);

%--------------------------------------------------------------------------
function [D, Do, h] = mixedMappings(G, W, D, Dw, h, hw, bc)
if ~isempty(W)
    nperf = cellfun(@numel, { W.cells });
else
    nperf = 0;
end
[fluxFacesR, fluxFacesW] = getFluxFaces(G, W, bc);
ncf = size(G.cellFaces, 1);
cellNo = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2)';
sgn    = 2*( G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo ) -1;
Do = blkdiag( spdiags(sgn, 0, ncf, ncf) * D, ...
              speye(sum(nperf), sum(nperf)) );
D  = blkdiag( D(:, fluxFacesR), Dw(:, fluxFacesW) );
h  = [h(fluxFacesR); hw(fluxFacesW)];

%--------------------------------------------------------------------------
function [fluxFacesR, fluxFacesW] = getFluxFaces(G, W, bc)
fluxFacesR = (prod(double(G.faces.neighbors), 2) == 0);
if ~isempty(bc)
    pf = strcmpi('pressure', bc.type');
    fluxFacesR ( double(bc.face(pf)) ) = false;
end
if nargout > 1
    if ~isempty(W)
        fluxFacesW = strcmpi('rate', {W.type}');
    else
        fluxFacesW = [];
    end
end

%--------------------------------------------------------------------------
function [resSol, wellSol] = packSol(G, W, resSol, wellSol, flux, p, lamN, bc)
fluxFacesR = getFluxFaces(G, W, bc);
facePressure = zeros(G.faces.num, 1);
if ~isempty(bc)
    pf = strcmpi('pressure', bc.type');
    presFaces = double(bc.face(pf));
    facePressure(presFaces) = bc.value(pf);
end
facePressure(fluxFacesR) = lamN(1:nnz(fluxFacesR));

resSol.facePressure = facePressure;
resSol.cellPressure = p;
resSol.faceFlux     = flux(1:G.faces.num);
resSol.cellFlux     = faceFlux2cellFlux(G, flux);

%wells
r_inx = G.faces.num;
p_inx = numel(fluxFacesR);
for k = 1:numel(W)
    if strcmpi(W(k).type, 'pressure')
        wellSol(k).pressure = W(k).val;
    else
        wellSol(k).pressure = lamN(p_inx+1);
        p_inx = p_inx+1;
    end
    nwc = numel(W(k).cells);
    wellSol(k).flux = - flux(r_inx + (1:nwc));
    r_inx = r_inx + nwc;
end
 
