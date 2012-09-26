function resSol = implicitTransport(resSol, wellSol, G, tf, ...
                                    rock, fluid, varargin)
%Implicit single point upwind transport solver for two-phase flow.
%
% SYNOPSIS:
%   resSol = implicitTransport(resSol, wellSol, G, tf, rock, fluid)
%   resSol = implicitTransport(resSol, wellSol, G, tf, rock, fluid, ...
%                              'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function implicitTransport solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf] by calling function initTransport and one of the
%   transport solvers twophaseUpwBE or twophaseUpwBEgrav depending on
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
%
%   OnlyGrav - Only consider transport caused by gravity, (ignore Darcy
%              flux from pressure solution).  Used for gravity splitting.
%              Default value: OnlyGrav = false.
%
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
% EXAMPLE:
%   See simple2phWellExample.m
%
% SEE ALSO:
%   initTransport, twophaseUpwBE, twophaseUpwBEGrav, explicitTransport.

%{
#COPYRIGHT#
%}

% $Id: implicitTransport.m 2964 2009-10-09 05:08:08Z jrn $

opt  = struct('verbose' , false , ...  % Emit progress reports?
              'nltol'   , 1.0e-6, ...  % Non-linear residual tolerance
              'lstrials', 20    , ...  % Max no of line search trials
              'maxnewt' , 25    , ...  % Max no. of NR iterations
              'tsref'   , 12    , ...  % Time step refinement
              'resred'  , 0.99  , ...  % Residual reduction factor
              'OnlyGrav', false , ...
              'wells'   , []    , ...
              'src'     , []    , ...
              'bc'      , []);

opt = merge_options(opt, varargin{:});

[g, flux, q, pv] = initTransport(G, resSol, wellSol,              ...
                                 rock, fluid, 'ComputeDt', false, ...
                                 'OnlyGrav', opt.OnlyGrav,        ...
                                 'wells', opt.wells, 'src', opt.src, ...
                                 'bc', opt.bc);

opt = rmfield(opt, {'OnlyGrav', 'src', 'bc', 'wells'});
arg = [fieldnames(opt), struct2cell(opt)] .';

g_vec = gravity();
if norm(g_vec(1 : size(G.nodes.coords,2))) > 0,
   [resSol, report] = twophaseUpwBEGrav(resSol, G, tf, q, flux, g, ...
                              pv, fluid, arg{:});
   dispReport(report);
else
   resSol = twophaseUpwBE(resSol, tf, q, flux, pv, fluid, arg{:}, ...
                          'LinSolve', @mldivide);
end

if any(any(isnan(resSol.s))),
   disp('Transport step failed')
end


function dispReport(r)
   if r.success,
      fprintf(['Transport step used %d iterations in %d sub steps ',...
               '(%d failed steps)\n'], r.iterations, r.sub_steps,   ...
                r.failed_steps);
   else
      fprintf(['Transport step FAILED. (%d iterations, %d sub steps, ',...
               '%d failed steps)\n'], r.iterations, r.sub_steps,   ...
                r.failed_steps);
   end   

