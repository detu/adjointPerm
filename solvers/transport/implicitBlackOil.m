function resSol = implicitBlackOil(G, resSol, wellSol, tf, pv, fluid, varargin)
%Implicit single point upwind transport solver for Black Oil flow.
%
% SYNOPSIS:
%   resSol = implicitBlackoil(G, resSol, wellSol, tf, pv, fluid)
%   resSol = implicitBlackoil(G, resSol, wellSol, tf, pv, fluid, ...
%                             'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function implicitBlackoil computes the solution to the BlackOil transport
%   equation
%
%        z_t + \div·(rho·vf(s))_x = q, with u = z / B(p) and s = u/sum(u),
%
%   using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time.  The system of nonlinear equations that must be
%   solved to move the solution from time=0 to time=t, is solved using a
%   Newton-Raphson algorithm with line search to increase robustness.
%
%   In the case of failing to compute a solution using only a single step
%   over the interval [0,t], an alternative strategy involving sub-steps
%   and step size control is employed.
%
% REQUIRED PARAMETERS:
%   G        - Grid data structure discretising the geometric model.
%
%   resSol   - Initial solution (at time=0).  The field z (density) must be
%              present with one scalar value for each phase in each cell in
%              the grid.
%
%   wellSol  - Well solution structure giving well fluxes and well
%              pressures at time=0.
%
%   tf       - Length of time step.
%
%   pv       - Reservoir pore volumes, one scalar value for each cell in the
%              grid.
%
%   fluid    - struct holding functions fluid.pvt and fluid.relperm.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   verbose  - Whether or not time integration progress should be reported
%              to the screen.
%              Default value: verbose = false.
%
%   nltol    - Absolute tolerance of iteration.  The numerical solution
%              must satisfy the condition
%
%                NORM(S-S0 + dt/porvol(out - in) - Q, INF) <= nltol
%
%              at all times in the interval [0,t].
%
%              Default value: nltol = 1.0e-3.
%
%   lstrials - Maximum number of trials in linesearch method.  Each new
%              trial corresponds to halving the step size along the search
%              direction.
%              Default value: lstrials = 20.
%
%   maxnewt  - Maximum number of inner iterations in Newton-Raphson method.
%              Default value: maxnewt = 25.
%
%   tsref    - Maximum time step refinement power.  The minimum time step
%              allowed is t / 2^tsref.
%              Default value: tsref = 12.
%
% RETURNS:
%   resSol   - Updated solution.  Fields z (density) and s (saturation) are
%              modified.
%
% SEE ALSO:
%   twophaseUpwFE, evalpvt.

%{
#COPYRIGHT#
%}

% $Id: implicitBlackOil.m 2338 2009-06-05 17:19:30Z bska $


   opt      = struct('verbose',  true,   ... % emit progress reports
                     'nltol',    1.0e-9, ... % non-linear residual tolerance
                     'lstrials', 20,     ... % max no of line search trials
                     'maxnewt',  25,     ... % max no. of NR iterations
                     'tsref',    12,     ... % time step refinement
                     'resred',   0.99,   ... % residual reduction factor
                     'bc',       [],     ...
                     'src',      [],     ...
                     'wells',    []);

   opt      = merge_options(opt, varargin{:});

   [F, Jac] = blackoilJacobian2(G, resSol, wellSol, pv, fluid, ...
                               opt.wells, opt.src, opt.bc);

   update   = @(zn, z0, dz, dt, err) linesearch(zn, dz, opt.resred * err, ...
                                                @(z) F(z, z0, dt),    ...
                                                opt.lstrials);
%   update2   = @(zn, z0, d, err) linesearch2(zn, d, opt.resred * err, ...
%                                              @(z, dt) F(z, z0, dt),    ...
%                                               opt.lstrials);

   z        = reshape(resSol.z', [], 1);
%  z        = ContinuationNewton(z, tf, F, Jac, JacT, update2, opt);
   z        = newtonRaphson(z, tf, F, Jac, update, opt);
   resSol.z = reshape(z, size(resSol.z'))';
   [u,u,u,u] = fluid.pvt(resSol.cellPressure, resSol.z);
   resSol.s = bsxfun(@rdivide, u, sum(u,2));
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%{
function [s, res, linfail, alph] = linesearch2(s, d, target_residual_reduction, F, number_of_trials)
   alph = 0;
   i    = 0;
   linfail = true;

   % Geometric line search s+ds, s+ds/2, s+ds/4, ...
   while linfail && (i < number_of_trials),
     da = pow2(d, alph);
     sn   = s + da(1:end-1);
     res  = F(sn, da(end));

     alph = alph - 1;
     i    = i + 1;
     linfail = ~(norm(res, inf) < target_residual_reduction);
   end

   %alph = pow2(alph + 1);      % Undo last (unneeded) scaling.
   s    = sn;
end
%}

function [s, res, linfail] = linesearch(s, ds, target_residual_reduction, F, number_of_trials)
   alph = 0;
   i    = 0;
   linfail = true;

   % Geometric line search s+ds, s+ds/2, s+ds/4, ...
   while linfail && (i < number_of_trials),
     sn   = s + pow2(ds, alph);
     res  = F(sn);

     alph = alph - 1;
     i    = i + 1;
     linfail = ~(norm(res, inf) < target_residual_reduction);
   end

   %alph = pow2(alph + 1);      % Undo last (unneeded) scaling.
   s    = sn;
end


function [s, res, linfail] = simplechop(s, ds, target_residual_reduction, F, number_of_trials)
  % Chop Newton update to ensure that saturations remain in [0,1] and that the
  % change is less than LIM.

  lim  = 0.2;
  ds   = min(S+ds, 1) - S;
  ds   = max(S+dx, 0) - S;
  ds   = max(min(ds, lim), -lim);
  s    = s+ds;
end
