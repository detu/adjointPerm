function resSol = twophaseUpwTheta(resSol, tf, q, gm, pv, fluid, varargin)
%Implicit single point upwind solver for two-phase flow, no gravity.
%
% SYNOPSIS:
%   resSol = twophaseUpwTheta(resSol, t, q, gm, porvol, fluid)
%   resSol = twophaseUpwTheta(resSol, t, q, gm, porvol, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwTheta solves the Buckley-Leverett transport equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a theta scheme
%   discretisation in time. The nonlinear system of equations that must be
%   solved to move the solution from time=0 to time=tf, are solved using a
%   Newton-Raphson algorithm with line search to increase robustness.
%
%   In the case of failing to compute a solution using only a single step
%   over the interval [0,tf], an alternative strategy involving sub-steps
%   and step size control is employed.
%
% REQUIRED PARAMETERS:
%   resSol   - Reservoir solution structure containing valid (water)
%              saturation resSol.s(:,1) with one value for each cell in the
%              reservoir model.
%
%   tf       - Length of time step, measured in units of seconds.
%
%   q        - Accumulated sources (typically contributions from wells).
%              One scalar value for each cell in the grid.  The source
%              rates are assumed to be measured in units of m^3/s.
%
%   gm       - Upwind/inflow matrix of fluxes into each cell, created
%              (e.g.) by function initTransport.  Specifically, gm(i,j) is
%              the (positive) flux from cell j to cell i.  The fluxes are
%              assumed to be measured in units of m^3/s.
%
%   porvol   - Reservoir pore volumes, measured in units of m^3, one scalar
%              value for each cell in the reservoir model.
%
%   fluid    - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   verbose  - Whether or not time integration progress should be reported
%              to the screen.
%              Default value: verbose = false.
%
%   theta    - Parameter for theta scheme, valid values are in [0,1].
%              Default value = 0.5.
%
%
%   nltol    - Absolute tolerance of iteration.  The numerical solution
%              must satisfy the condition
%
%                NORM(S-S0 + dt/porvol(out - in) - Q, INF) <= nltol
%
%              at all times in the interval [0,t].
%
%              Default value: nltol = 1.0e-6.
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
%   resSol   - Reservoir solution with updated resSol.s.
%
% SEE ALSO:
%   twophaseUpwBEGrav, initTransport, implicitTransport, explicitTransport.

%{
#COPYRIGHT#
%}

% $Id: twophaseUpwTheta.m 2044 2009-04-20 17:18:53Z bska $

[F, Jac, linsrch, opt, dt] = setup(tf, q, gm, pv, fluid, varargin{:});
             
% use dt? 
resSol = newtonRaphson2ph(resSol, tf, F, Jac, linsrch, opt);

end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [F, Jac, linsrch, prm, dt] = setup(tf, q, gm, pv, fluid, varargin)

nc = numel(pv);

if ~((ndims(gm) == 2) && all(size(gm) == nc)),
    error('twophaseUpwTheta:SETUP:WRONG_SIZE',           ...
         ['Unexpected size of inflow flux matrix.\n', ...
          'Expected [', int2str([nc, nc]), '],\n',    ...
          'but got  [', int2str(size(gm)), ']']);
end

prm  = struct(  'verbose',  false,  ...  % emit progress reports
                'nltol',    1.0e-6, ...  % non-linear residual tolerance
                'lstrials', 20,     ...  % max no of line search trials
                'maxnewt',  25,     ...  % max no. of NR iterations
                'tsref',    12,     ...  % time step refinement
                'resred',   0.99,   ...  % residual reduction factor
                'theta',    0.5,    ...  % default: Crank-Nicolson Scheme 
                'dt',       tf,     ...
                'LinSolve', @mldivide);
             
prm = merge_options(prm, varargin{:});

theta   = prm.theta;
dt = prm.dt;
[im, jm, valm] = find(gm);

q  = q (:);
pv = pv(:);

% System F(s) = 0 of non-linear equations (and the system's Jacobian
% matrix) defining saturation equilibrium at a given time step.
%
% Here, we use a Buckley-Leverett model with
%
%            s²/µw              s²                µw
%    f(s) = ---------------- = ------------ ,  mr=---
%            s²/µw+(1-s)²/µo    s²+mr*(1-s)²      µo
%
%
%            2 mr*s(1-s)
%   df(s) = ---------------
%           (s²+mr*(1-s)²)²
%
%
% With the matrix v, we can write an upwind theta-scheme
% discretisation of the Buckley-Leverett model as
%
%    s^n-s^{n-1} + dt*v*f(s) = MAX(q,0),
%
% where
%
%    v(i,j) = - flux into cell i from cell j,
%    v(i,i) = sum of flux out of cell i -MIN(q(i),0).
%
% The target function is
%
%    F(s) = s-s^{n-1} + dt*v*fs(s) - max(q, 0)
%
% and the Jacobian is
%
%   dF(s) =  I + dt*v*diag(df(s)).
%
f    = @(kr) (  kr(:,1)./(kr(:,1) + (fluid.mu(1)./fluid.mu(2))*kr(:,2)) );
df_w = @(kr, dkr, Lt) ((dkr(:,1)./fluid.mu(1))./Lt - ...
                (kr(:,1)./fluid.mu(1)).* ...
                sum(bsxfun(@rdivide, dkr, fluid.mu),2)./(Lt.^2));

% v <- (gm - colsum(gm) + MIN(q,0))) / pore_volume
v = spdiags(1./pv, 0, nc, nc) * sparse( [  im;    jm;  (1 : nc)'], ...
                                        [  jm;    jm;  (1 : nc)'], ...
                                        [valm; -valm; -min(q,0) ], nc, nc);
                                     
F   = @(sol, s0, dt, expl) (sol.s(:,1) - s0 + dt.* ...
                     (v*(theta*f(fluid.kr(sol)) + ...
                     (1-theta)*f(fluid.kr(struct('s', s0)))) - max(q./pv,0)));
Jac = @(sol, dt) ( speye(nc) + dt*theta.*v*sparse(1:nc, 1:nc, ...
                  df_w( fluid.kr(sol), fluid.dkr(sol), fluid.Lt(sol) ), ...
                  nc, nc, nc) );
              

% Bind prm.resred * err, @(sat) F(sat, s0, dt) and prm.lstrials
% to function call to reduce clutter in main code.
linsrch = @(sol, s0, ds, dt, err) linesearch(sol, ds, prm.resred * err, ...
                                             @(sol) F(sol, s0, dt),     ...
                                             prm.lstrials);
                                          
prm = rmfield(prm, {'lstrials', 'resred'});
end

%--------------------------------------------------------------------------

function [sol, res, alph, fail] = linesearch(sol, ds, target, F, ni)
%
% Basic idea: search for a step size 'alpha' in direction 'ds', subject to
% the restriction that alpha be in [0,1], which ensures that the objective
% function 'F' decreases.  That is: F(s + alpha*ds) < F(s).
%
% In the current implementation, alpha is reduced in a geometric sequence.
% A more sophisticated approach would ensure a certain minimum reduction as
% well.
minSat = 0;     % Minimum (water) saturation
maxSat = 1;     % Maximum (water) saturation
capSat = @(sat) min(max(minSat, sat), maxSat);

alph = 0;
i    = 0;
fail = true;

% Geometric line search: seems pretty robust
while fail && (i < ni),
   sn.s   = capSat(sol.s(:,1) + pow2(ds, alph));
   res  = F(sn);

   alph = alph - 1;
   i    = i + 1;
   fail = ~(norm(res, inf) < target);
end

alph       = pow2(alph + 1);      % Undo last (unneeded) scaling.
sol.s(:,1) = sn.s;
end
