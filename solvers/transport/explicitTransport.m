function resSol = explicitTransport(resSol, wellSol, G, tf, ...
                                    rock, fluid, varargin)
%Explicit single point upwind transpot solver for two-phase flow.
%
% SYNOPSIS:
%   resSol = explicitTransport(resSol, wellSol, G, tf, rock, fluid)
%   resSol = explicitTransport(resSol, wellSol, G, tf, rock, fluid,
%                              'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function explicitTransport solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf] by calling function initTransport and one of the
%   transport solvers twophaseUpwFE or twophaseUpwFEgrav depending on
%   presence or absence of effects of gravity in the model.  See these
%   functions for additional details.
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
%   G       - Grid data structure discretising the reservoir model.
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
%   verbose   - Whether or not to emit time integration progress.
%               Default value: verbose = false.
%
%   wells     - Well structure as defined by functions 'addWell' and
%               'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%               which is interpreted as a model without any wells.
%
%   bc        - Boundary condtion structure as defined by function
%               'addBC'. This structure accounts for all external boundary
%               contributions to the reservoir flow.
%               Default value: bc = [] meaning all external no-flow
%               (homogeneous Neumann) conditions.
%
%   src       - Explicit source contributions as defined by function
%               'addSource'.  Default value: src = [] meaning no explicit
%               sources exist in the model.
%
%   OnlyGrav  - Only consider transport caused by gravity, (ignore darcy
%               flux from pressure solution). Used for gravity splitting.
%               Default value: OnlyGrav = false.
%
%   ComputeDt - Whether or not to compute timestep honouring CFL
%               condition. Default value: ComputeDt = true.
%
%   max_dt    - Maximal timestep constraint (measured in seconds).
%               Default value: max_dt = 100 * day() (100 days).
%
%   dt_factor - Factor to multiply dt by.
%               Default value: dt_factor = 0.5.
%
%   dt        - User specified internal timestep (measured in seconds).
%               Ignored if ComputeDt = true and dt > computed time step.
%               Default value: dt = tf.
%
%               NB: The explicit scheme is only stable provided that dt
%               honours a CFL time step restriction.
%
%   satWarn   - Tolerance level for saturation warning.
%               Default value: satWarn = sqrt(eps).
%
% RETURNS:
%   resSol - Reservoir solution with updated saturation, resSol.s.
%
% EXAMPLE:
%   twophaseTransport
%
% SEE ALSO:
%   initTransport, twophaseUpwFE, twophaseUpwFEGrav, implicitTransport.

%{
#COPYRIGHT#
%}

% $Id: explicitTransport.m 2200 2009-05-15 14:34:54Z hnil $

opt = struct('verbose'  , false      , ...  % Emit progress reports?
             'OnlyGrav' , false      , ...
             'ComputeDt', true       , ...
             'max_dt'   , 100 * day(), ...  % 100 days (in seconds)
             'dt_factor', 0.5        , ...
             'wells'    , []         , ...
             'src'      , []         , ...
             'bc'       , []         , ...
             'dt'       , tf         , ...
             'satWarn'  , sqrt(eps));

opt = merge_options(opt, varargin{:});

% Verify that function will give valid internal timestep dt > 0.
dt = min(opt.dt, opt.max_dt);
assert (opt.ComputeDt || dt > 0);

verbose = opt.verbose;
satWarn = opt.satWarn;
opt     = rmfield(opt, {'dt','satWarn'});
arg = [fieldnames(opt), struct2cell(opt)].';

[g, flux, q, pv, est_dt] = initTransport(G, resSol, wellSol, ...
                                         rock, fluid, arg{:});

% Set internal timestep.
if opt.ComputeDt,
   dt = min(dt, est_dt);
end

g_vec = gravity();
if norm(g_vec(1 : size(G.nodes.coords,2))) > 0,
   resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, g, pv, fluid, ...
                              'dt', dt, 'verbose', verbose, ...
                              'satWarn', satWarn);
else
   resSol = twophaseUpwFE(resSol, tf, q, flux, pv, fluid, 'dt', dt, ...
                          'verbose', verbose, 'satWarn', satWarn);
end
