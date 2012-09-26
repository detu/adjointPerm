function [resSol,report] = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, ...
                                    pv, fluid, varargin)
%Implicit single point upwind solver for two-phase flow, including gravity.
%
% SYNOPSIS:
%   resSol = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, pv, fluid)
%   resSol = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, pv, fluid, ...
%                              'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwBEGrav solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time.  The nonlinear system of equations that must be
%   solved to move the solution from time=0 to time=tf, are solved using a
%   Newton-Raphson algorithm with line search to increase robustness.
%
%   In the case of failing to compute a solution using only a single step
%   over the interval [0,tf], an alternative strategy involving sub-steps
%   and step size control is employed.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   q       - Accumulated sources (typically contributions from wells).
%             One scalar value for each cell in the grid.  The source rates
%             are assumed to be measured in units of m^3/s.
%
%   flux    - Inflow matrix of fluxes into each cell, created
%             by function initTransport.  Size
%             G.cells.num-by-(G.faces.num-boundary faces).  The fluxes are
%             assumed to be measured in units of m^3/s.
%
%   grav    - Matrix with gravity contribution for each face, created by
%             function initTransport.  Same size as flux.
%
%   pv      - Reservoir pore volumes, measured in units of m^3.  One scalar
%             value for each cell in the model.
%
%   fluid   - Data structure describing the fluids in the problem.
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
%              allowed is tf / 2^tsref.
%              Default value: tsref = 12.
%
% RETURNS:
%   resSol   - Reservoir solution with updated resSol.s.
%
% SEE ALSO:
%   twophaseUpwBE, initTransport, implicitTransport, explicitTransport.

% TODO:
%   - implement gravity effects for pressure boundary and wells

%{
#COPYRIGHT#
%}

% $Id: twophaseUpwBEGrav.m 2964 2009-10-09 05:08:08Z jrn $

assert (size(resSol.s,2)<3 || all(resSol.s(:,3)==0));

[F, Jac, linsrch, verbose, ...
 nltol, maxnwt, mints] = setup(tf, q, flux, grav, pv, fluid, varargin{:});

[v_darcy, neighbors, ...
 normals, up_inx, grav_vec] = initFaceMob(G, resSol, flux, grav);

findMobM = @(mob) findFaceMobMat(up_inx, mob, grav_vec, v_darcy, ...
                                 neighbors, normals, fluid);
                              
report = struct('success',      true,...
                'iterations',   0, ...
                'sub_steps',    0, ...
                'failed_steps', 0);
                              
                              
sz = size(resSol.s(:,1));
sn = resSol;

dispif(verbose, '\n\n');
dispif(verbose, [repmat('-', [1, 70]), '\n']);
dispif(verbose, '  Time interval           iter   relax    residual \n');
dispif(verbose, [repmat('-', [1, 70]), '\n']);

[t, dt, dtprev, count] = deal(0.0, tf, -tf, 0);

%--------------------------------------------------------------------------
% Main time loop (ideally, just a single step) ----------------------------
%
while t < tf && dt > mints,
   dt = min(dt, tf - t);

   % Update faceMob index matrices:
   % NB: A_w, and A_o could be updated inside the inner NR algorithm
   % (before computing the Jacobian matrix) to increase accuracy, but doing
   % so would decrease the efficiency greatly.
   if t == 0 || nnz(grav)
      mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
      [A_w, A_o] = findMobM(mob);
   end

   %-----------------------------------------------------------------------
   % Outer controlling loop (sub step size &c) ----------------------------
   %
   redo_newton = true;
   while redo_newton,
      dispif(verbose, '[%f, %f]:', t / tf, (t + dt) / tf); nc = 0;

      s0 = resSol.s(:,1);
      sn = resSol; sn.s(:,1) = 0.5;

      %--------------------------------------------------------------------
      % Initialise NR algorithm -------------------------------------------
      %
      res     = F(sn, s0, dt, A_w, A_o);

      err     = norm(res, inf);
      nwtfail = err > nltol;

      linfail = false;
      it      = 0;
      %-----------------------------------------------------------------
      % Execute inner NR algorithm -------------------------------------
      %
      while nwtfail && ~linfail && it < maxnwt,

         ds = - Jac(sn, dt, A_w, A_o) \ res;

         [sn, res, alph, linfail] = linsrch(sn, s0, ds, dt, err, A_w, A_o);

         it      = it + 1;
         err     = norm(res, inf);
         nwtfail = err > nltol;

         dispif(verbose, repmat('\b', [1, nc]));
         nc = dispif(verbose, '\t%4d\t %2.2f\t %5.5e', it, alph, err);
         report.iterations = report.iterations +1;
      end

      %-----------------------------------------------------------------
      % Determine if NR succeeded or not -------------------------------
      %
      count = count - 1;
      if nwtfail,
         dispif(verbose, '\tReducing step.');

         if count > 0 && dtprev > 0,
            dt = dtprev;
         else
            dt = dt / (1.5 + 0.5);
         end
         count = 5;
         report.failed_steps = report.failed_steps + 1;
      else
         redo_newton = false;
         t = t + dt; dtprev = -dt;
         if it == 0,
            dispif(verbose, '\t%4d\t %s\t %5.5e', it, '-', err);
            dispif(verbose, '\t NB: err <= ntol.');
         end

         if it <= 5  && count < 0 && t < tf, % Arbitrary threshold.
            dispif(verbose, '\tIncreasing step.');

            dtprev = dt;
            dt     = min((1 + 0.5) * dt, tf - t);
            count  = 5;
         end
      end
      dispif(verbose, '\n');
   end

   % Time step [t, t+dt] was successful
   report.sub_steps = report.sub_steps + 1;
   resSol = sn;
end

if ~(t < tf),
   dispif(verbose, 'We''re done (%f)\n', t);
else
   dispif(verbose, 'Unable to integrate to final time, t=%f\n', tf);
   resSol.s(:,1) = NaN(sz);
   report.success = false;
end

if size(resSol.s,2) > 1,
   % Update oil saturation:
   resSol.s(:,2) = 1 - resSol.s(:,1);
end
end
%-----------------------------------------------------------------------
% Private helpers follow.
%-----------------------------------------------------------------------

function [F, Jac, linsrch, verbose, nltol, maxnwt, mints] = ...
      setup(tf, q, flux, grav, pv, fluid, varargin)

nc = numel(pv);
nf = size(flux,2);

prm  = struct(  'verbose',  false,  ...  % emit progress reports
                'nltol',    1.0e-6, ...  % non-linear residual tolerance
                'lstrials', 20,     ...  % max no of line search trials
                'maxnewt',  25,     ...  % max no. of NR iterations
                'tsref',    12,     ...  % time step refinement
                'resred',   0.99);       % residual reduction factor

prm = merge_options(prm, varargin{:});

verbose = prm.verbose;
nltol   = prm.nltol;
maxnwt  = prm.maxnewt;
mints   = pow2(tf, -prm.tsref);        % minimum time step size

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
% With the matrices flux, grav, A_w and A_o, we can write an upwind
% backward Euler discretisation of the Buckley-Leverett model as
%
%     s^(n+1) = s^n - (dt./pv)*(H(s^n+1) - max(q,0) - min(q,0)*f(s^n+1))
%
% where H(s) = (flux + grav*diag(A_o*lam_o(s))
%                 *A_w*lam_w(s)./(A_w*lam_w(s)+A_o*lam_o(s)),
%
% lam_l is the mobility for phase l, f is Buckely-Leverett fractional
% flow function, while A_o and A_w are index matrices that determine
% the upstream mobility.
%
% The target function is
%
%    F(s) = s-s^{n-1} + (dt./pv)*(H(s) - max(q, 0)- min(q,0)*f(s))
%
% and the Jacobian is
%
%   dF(s) =  I + (dt./pv)(H'(s)  - min(q,0)*df(s))

df_w = @(kr, dkr, Lt) ((dkr(:,1)./fluid.mu(1))./Lt -  ...
                        (kr(:,1)./fluid.mu(1)).*      ...
                        sum(bsxfun(@rdivide, dkr, fluid.mu),2)./(Lt.^2));

diaglam_o = @(A_o, kr)(sparse(1:nf,1:nf, A_o*kr(:,2)./fluid.mu(2),nf, nf));

   function [kr, face_den, f_face] = compVec(resSol, A_w, A_o)
      kr       = fluid.kr(resSol);
      face_den = (A_w*kr(:,1) + (fluid.mu(1)/fluid.mu(2)).*A_o*kr(:,2));
      f_face   = (A_w*kr(:,1))./face_den;
   end

   function F = compF(resSol, s0, dt, A_w, A_o)
      [kr, face_den, f_face] = compVec(resSol, A_w, A_o);

      f = kr(:,1)./(kr(:,1)+(fluid.mu(1)/fluid.mu(2)).*kr(:,2));
      F = resSol.s(:,1) - s0 + dt.*(1./pv).*( ...
          (flux +  grav * diaglam_o(A_o, kr))*f_face...
          - max(q,0) - min(q,0).*f);
   end

   function J = compJacGrav(resSol, dt, A_w, A_o)
      dkr = fluid.dkr(resSol);
      Lt  = fluid.Lt(resSol);

      C = (fluid.mu(1)/fluid.mu(2));
      [kr, face_den, f_face] = compVec(resSol, A_w, A_o);
      dw      = A_w *sparse(1:nc, 1:nc, dkr(:,1), nc, nc);
      do      = A_o* sparse(1:nc, 1:nc, dkr(:,2), nc, nc);
      df_face = sparse(1:nf, 1:nf, 1./face_den)*dw ...
                - sparse(1:nf, 1:nf ,f_face./face_den)*(dw+do.*C);

      J = speye(nc) + sparse(1:nc, 1:nc, dt*(1./pv))*     ...
          ((flux + grav*diaglam_o(A_o, kr))*df_face       ...
          + grav*sparse(1:nf,1:nf,f_face)*do./fluid.mu(2) ...
          - sparse(1:nc,1:nc, min(q, 0))*                 ...
          sparse(1:nc, 1:nc, df_w(kr, dkr, Lt)));
   end

   function J = compJacNoGrav(resSol, dt, A_w)
      kr  = fluid.kr(resSol);
      dkr = fluid.dkr(resSol);
      Lt  = fluid.Lt(resSol);

      J  = speye(nc) + sparse(1:nc, 1:nc,dt*(1./pv), nc, nc) ...
         *(flux*A_w - sparse(1:nc,1:nc, min(q, 0)))...
         *sparse(1:nc, 1:nc, df_w(kr, dkr, Lt), nc, nc, nc);
   end


   % Return function handles
   F = @compF;
   if   nnz(grav), Jac = @compJacGrav;
   else            Jac = @(varargin) compJacNoGrav(varargin{1:3});
   end

   % Bind prm.resred * err, @(sat) F(sat, s0, dt) and prm.lstrials
   % to function call to reduce clutter in main code.
   linsrch = @(resSol, s0, ds, dt, err, A_w, A_o)        ...
               linesearch(resSol, ds, prm.resred * err,  ...
               @(resSol) F(resSol, s0, dt, A_w, A_o), prm.lstrials);
end

%--------------------------------------------------------------------------

function [resSol, res, alph, fail] = linesearch(resSol, ds, target, F, ni)
%
% Basic idea: search for a step size 'alpha' in direction 'ds', subject to
% the restriction that alpha be in [0,1], which ensures that the objective
% function 'F' decreases.  That is: F(s + alpha*ds) < F(s).
%
% In the current implementation, alpha is reduced in a geometric sequence.
% A more sophisticated approach would ensure a certain minimum reduction as
% well.
%
   minSat = 0;     % Minimum (water) saturation
   maxSat = 1;     % Maximum (water) saturation
   capSat = @(sat) min(max(minSat, sat), maxSat);

   alph = 0;
   i    = 0;
   fail = true;

   % Geometric line search: seems pretty robust
   while fail && (i < ni),
      sn.s = capSat(resSol.s(:,1) + pow2(ds, alph));
      res  = F(sn);

      alph = alph - 1;
      i    = i + 1;
      fail = ~(norm(res, inf) < target);
   end

   alph = pow2(alph + 1);      % Undo last (unneeded) scaling.
   resSol.s(:,1) = sn.s;
end
