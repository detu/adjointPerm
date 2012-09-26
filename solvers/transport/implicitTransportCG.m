function resSol = implicitTransportCG(resSol, wellSol, G, CG, tf, ...
                                      rock, fluid, varargin)
% Coarse grid implicit single point upwind transport for two-phase flow
%
% SYNOPSIS:
%   resSol = implicitTransportCG(resSol, wellSol, G, CG, tf, rock, fluid)
%   resSol = implicitTransportCG(resSol, wellSol, G, CG, tf, rock, fluid, ...
%                                'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function implicitTransportCG is the coarse grid version of function
%   implicitTransport, see function for details. The transport problem is
%   solved on a coarse saturation grid which need not be the same grid as
%   used for the pressure solution. The coarse saturation is then projected
%   onto the fine (pressure) grid, and the fine scale reservoir solution
%   structure is returned.  
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid (water)
%             saturation resSol.s(:,1) with one value for each cell
%             in the grid.  Pressures are assumed to be measured in units
%             of Pascal while fluxes are assumed to be measured in units of
%             m^3/s.
%
%   wellSol - Well solution structure.  Pressures are assumed to be
%             measured in units of Pascal while fluxes are assumed to be
%             measured in units of m^3/s.
%
%   G       - Grid data structure discretising the reservoir model
%             (pressure grid).
%
%   CG      - Coarse grid (saturation grid).
%
%   tf      - End point of time integration interval (i.e., final time).
%             Measured in units of seconds.
%
%   rock    - Rock data structure.  Must contain the field 'rock.poro',
%             and in the presence of gravity, valid permeabilities measured
%             in units of m^2 in field 'rock.perm'.
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   verbose  - Whether or not time integration progress should be
%              reported to the screen. Default value: verbose = false.
%
%   wells    - Well structure as defined by functions 'addWell' and
%              'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%              which is interpreted as a model without any wells.
%
%   bc       - Boundary condtion structure as defined by function
%              'addBC'. This structure accounts for all external boundary
%              contributions to the reservoir flow.
%              Default value: bc = [] meaning all external no-flow
%              (homogeneous Neumann) conditions.
%
%   src      - Explicit source contributions as defined by function
%              'addSource'. Default value: src = [] meaning no explicit
%              sources exist in the model.
%   nltol    - Absolute tolerance of iteration.  The numerical solution
%              must satisfy the condition
%
%                 NORM(S-S0 + dt/porvol(out - in) - Q, INF) <= nltol
%
%              at all times in the interval [0,tf].
%              Default value: nltol = 1.0e-6.
%
%   lstrials - Maximum number of trials in linesearch method.  Each new
%              trial corresponds to halving the step size along the
%              search direction. Default value: lstrials = 20.
%
%   maxnewt  - Maximum number of inner iterations in Newton-Raphson method.
%              Default value: maxnewt = 25.
%
%   tsref    - Maximum time step refinement power.  The minimum time step
%              allowed is tf / 2^tsref.
%              Default value: tsref = 12.
%
% RETURNS:
%   resSol - Reservoir solution with updated saturation, resSol.s.
%
%
% SEE ALSO:
%   initTransport, twophaseUpwBE, twophaseUpwBEGrav, implicitTransport,
%   explicitTransport.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-14 12:13:58 +0200 (on, 14 okt 2009) $
% $Revision: 3003 $

opt  = struct('verbose' , false , ...  % Emit progress reports?
              'nltol'   , 1.0e-6, ...  % Non-linear residual tolerance
              'lstrials', 20    , ...  % Max no of line search trials
              'maxnewt' , 25    , ...  % Max no. of NR iterations
              'tsref'   , 12    , ...  % Time step refinement
              'resred'  , 0.99  , ...  % Residual reduction factor
              'wells'   , []    , ...
              'src'     , []    , ...
              'bc'      , []);

opt = merge_options(opt, varargin{:});

[g, flux, q, pv] = initTransport(G, resSol, wellSol,              ...
                                 rock, fluid, 'ComputeDt', false, ...
                                 'OnlyGrav', false,        ...
                                 'wells', opt.wells, 'src', opt.src, ...
                                 'bc', opt.bc);
                              
% ----------------------------------------------------------------------- 
% Coarse version 
% ----------------------------------------------------------------------- 

% Projection matrix from fine pressure grid to coarse saturation grid
Rc = CG.cells.subCells';

flux_CG = Rc*flux*Rc';
pv_CG   = Rc*pv;
q_CG    = Rc*q;

% Make coarse "solution structure" as input to twophaseUpwBE
resSol_CG = [];

% Add coarse saturation field
if isfield(resSol, 's_c')
   assert( numel(resSol.s_c) == CG.cells.num )
   resSol_CG.s = resSol.s_c;
else % convert saturation values from fine to coarse grid
   resSol_CG.s = fine2coarse(resSol.s, CG, G, rock); 
end
% ----------------------------------------------------------------------- 

opt = rmfield(opt, {'src', 'bc', 'wells'});
arg = [fieldnames(opt), struct2cell(opt)] .';

resSol_CG = twophaseUpwBE(resSol_CG, tf, q_CG, flux_CG, pv_CG, fluid, ...
                          arg{:}, 'LinSolve', @mldivide);

resSol.s_c = resSol_CG.s;
 
% Project coarse saturation onto fine grid 
resSol.s = coarse2fine(resSol_CG.s, CG);

if any(any(isnan(resSol.s))),
   disp('Transport step failed')
end

end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function sat_data = coarse2fine(sat_data_cg, CG)
%
% SYNOPSIS:
%   sat_data = coarse2fine(sat_data_cg, CG);
%
% DESCRIPTION:
%   Projects coarse saturation values to fine grid representation.
%
% PARAMETERS:
%
%   sat_data_cg - Vector of block saturation (saturation in coarse grid
%                 cells).
%   CG  - Coarse grid data structure. Assumes existing field:
%         CG.cells.subCells. Must be fixed when this field does not longer
%         exits.
%
% RETURNS:
%   sat_data - A vector of saturation values on fine grid.

sat_data = CG.cells.subCells * sat_data_cg;
end

function sat_data_cg = fine2coarse(sat_data, CG, G, rock)
%
% SYNOPSIS:
%   sat_data_cg = fine2coarse(sat_data, CG, rock, G);
%
% DESCRIPTION:
%    Computes coarse saturation values from fine grid representation.
%
%    First computes sum of volume of fluid in the block by multiplying the
%    fine saturation values with corresponding cell pore volmue, then
%    divide the total volume of fluid by the sum of pore volume in block.
%
% PARAMETERS:
%   sat_data - Vector of cell saturation. 
%   CG       - Coarse grid data structure. Assumes existing field:
%              CG.cells.subCells. Must be fixed when this field does not
%              longer exits.
%   G        - Grid data structure. rock - Rock structure.
%
% RETURNS:
%   sat_data_cg - A vector of saturation values on coarse grid blocks.

pv = poreVolume(G,rock);
sat_data_cg = zeros(CG.cells.num,1);

for i=1:CG.cells.num,
    cells = CG.cells.subCells(:,i);
    cells = find(cells); % dette kan kanskje gjøres bedre?
    sat_data_cg(i) = sum(sat_data(cells).*pv(cells))/sum(pv(cells));
end
end




