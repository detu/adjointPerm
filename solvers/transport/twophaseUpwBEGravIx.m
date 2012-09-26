function resSol = twophaseUpwBEGravIx(resSol, G, tf, q, flux, grav, ...
                                    pv, fluid, rock, varargin)
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

% $Date: 2009-06-08 13:05:26 +0200 (ma, 08 jun 2009) $
% $Revision: 2347 $

assert (size(resSol.s,2)<3 || all(resSol.s(:,3)==0));

[F, Jac, linsrch, verbose, ...
 nltol, maxnwt, mints] = setup(G,tf, q, flux, pv, fluid, rock, varargin{:});

[v_darcy, neighbors, ...
 normals, up_inx, grav_vec] = initFaceMob(G, resSol, flux, grav);

findMobIx = @(mob) findFaceMobIx(up_inx, mob, grav_vec, v_darcy, ...
                                 neighbors, normals, fluid);
                              
                              
%    function res = residual (sn, s0, dt)
%       mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
%       [iw, io] = findMobM(mob);
%       res     = F(sn, s0, dt, iw, io);
%    end
% 
%    function jac = jacobian (sn, dt)
%       mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
%       [iw, io] = findMobM(mob);
%       jac = Jac(sn, dt, iw, io);
%    end


sz = size(resSol.s(:,1));
sn = resSol;

%z = newtonRaphson(z, tf, F, Jac, update, opt)


if 1
   printif(verbose, '\n\n');
   printif(verbose, [repmat('-', [1, 70]), '\n']);
   printif(verbose, '  Time interval           iter   relax    residual \n');
   printif(verbose, [repmat('-', [1, 70]), '\n']);

   [t, dt, dtprev, count] = deal(0.0, tf, -tf, 0);

   %--------------------------------------------------------------------------
   % Main time loop (ideally, just a single step) ----------------------------
   %
   if verbose, tic, end;
   while t < tf && dt > mints,
      dt = min(dt, tf - t);

      %%{
   % Update faceMob index matrices:
   % NB: iw, and io could be updated inside the inner NR algorithm
   % (before computing the Jacobi) to increase accuracy, but doing so would
   % decrease the efficiency greatly.
   if t == 0 || nnz(grav)
      mob = bsxfun(@rdivide, fluid.kr(resSol), fluid.mu);
      [iw, io] = findMobIx(mob);
   end
      %}
      %-----------------------------------------------------------------------
      % Outer controlling loop (sub step size &c) ----------------------------
      %
      redo_newton = true;
      while redo_newton,
         printif(verbose, '[%f, %f]:', t/tf, (t + dt)/tf); nc = 0;

         s0 = resSol.s(:,1);
         sn = resSol; sn.s(:,1)=0.5;

         %--------------------------------------------------------------------
         % Initialise NR algorithm -------------------------------------------
         %
         %mob = bsxfun(@rdivide, fluid.kr(sn), fluid.mu);
         %[iw, io] = findMobIx(mob);
         res     = F(sn, s0, dt, iw, io);

         err     = norm(res, inf);
         nwtfail = err > nltol;

         linfail = false;
         it      = 0;
         %-----------------------------------------------------------------
         % Execute inner NR algorithm -------------------------------------
         %
         while nwtfail && ~linfail && it < maxnwt,

            %mob = bsxfun(@rdivide, fluid.kr(sn), fluid.mu);
            %[iw, io] = findMobIx(mob);
            ds = - Jac(sn, dt, iw, io) \ res;

            [sn, res, alph, linfail] = linsrch(sn, s0, ds, dt, err, iw, io);

            it      = it + 1;
            err     = norm(res, inf);
            nwtfail = err > nltol;

            printif(verbose, repmat('\b', [1, nc]));
            nc = printif(verbose, '\t%4d\t %2.2f\t %5.5e', it, alph, err);

         end

         %-----------------------------------------------------------------
         % Determine if NR succeeded or not -------------------------------
         %
         count = count - 1;
         if nwtfail,
            printif(verbose, '\tReducing step.');

            if count>0 && dtprev>0
               dt = dtprev;
            else
               dt    = dt / (1.5); %+rand(1));
            end
            count = 5;
         else
            redo_newton = false;
            t = t + dt; dtprev=-dt;
            if it == 0
               printif(verbose, '\t%4d\t %s\t %5.5e', it, '-', err);
               printif(verbose, '\t NB: err <= ntol.');
            end

            if it <= 5  && count < 0 && t ~= tf, % arbitrary threshold
               printif(verbose, '\tIncreasing step.');

               dtprev = dt;
               dt = min((1+rand(1)) * dt, tf - t);
               count = 5;
            end

         end
         printif(verbose, '\n');
      end

      % Time step [t, t+dt] was successful
      resSol = sn;

   end
   tocif(verbose)
end

if ~(t < tf),
   printif(verbose, 'We''re done (%f)\n', t);
else
   printif(verbose, 'Unable to integrate to final time, t=%f\n', tf);
   resSol.s(:,1) = NaN(sz);
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
   setup(G, tf, q, flux, pv, fluid, rock, varargin)

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
%    F(S) = S-S^{n-1} + (dt./pv)*(H(S) - max(q, 0)- min(q,0)*f(S))
%
% and the Jacobian is
%
%   dF(S) =  I + (dt./pv)(H'(S)  - min(q,0)*df(S))

   
   function F = compF(resSol, s0, dt, iw, io)
      i   = all(G.faces.neighbors~=0, 2);
      ic  = reshape(double(G.faces.neighbors(i,:)), [], 1);
      sgn = [ones(sum(i), 1); -ones(sum(i),1)];
      io  = reshape(repmat(double(io), [1, 2]), [], 1);
      iw  = reshape(repmat(double(iw), [1, 2]), [], 1);
      [g, f] = getFlux(G, rock, fluid, resSol, false);
      if ~any(g), assert(iw == io);end
      
      mob = fluid.mob(resSol);
      f_w_cell = mob(:,1) ./ sum(mob,2);
      
      mob_face = [mob(iw, 1) mob(io,2)];
      f_w_face = mob_face(:,1)./sum(mob_face,2);
      
      F = resSol.s(:,1) - s0 -dt.*(1./pv).*(max(q,0) + min(q,0).*f_w_cell);
      F = F + dt./pv.*accumarray(ic, sgn.*f.*f_w_face + ...
                                     sgn.*g.*(mob_face(:,2).*f_w_face));
      
   end

   function J = compJacGrav(resSol, dt, iw, io)
      i   = all(G.faces.neighbors~=0, 2);
      ic  = reshape(double(G.faces.neighbors(i,:)), [], 1);
      sgn = [ones(sum(i), 1); -ones(sum(i),1)];
      ic1 = (1:G.cells.num)';
      io  = reshape(repmat(double(io), [1, 2]), [], 1);
      iw  = reshape(repmat(double(iw), [1, 2]), [], 1);
      [g, f] = getFlux(G, rock, fluid, resSol, false);
      if ~any(g), assert(iw == io);end
            
      % [lam_w lam_o]
      mob = fluid.mob(resSol);
      % [(d/ds) lam_w (d/ds) lam_o]
      dkr_mu = bsxfun(@rdivide, fluid.dkr(resSol),fluid.mu);

      % Face flux
      mob_face = [mob(iw, 1) mob(io,2)];
      f_face   = bsxfun(@rdivide, mob_face, sum(mob_face,2));
      dw_fo    = f_face(:,2).*dkr_mu(iw,1);
      do_fw    = f_face(:,1).*dkr_mu(io,2);

      % Well flux out is computed in cell.
      %dfw = df_w(dkr_mu(ic1, :), mob(ic1, :));
      ff   = bsxfun(@rdivide, mob, sum(mob,2));
      df  = (ff(:,2).*dkr_mu(:,1) - ff(:,1).*dkr_mu(:,2))./sum(mob,2);
      
      
      zd = 1-dt./pv(ic1).*(min(q(ic1), 0).*df);
      zw = dt./pv(ic).*( sgn.*(f + g.*mob(io,2)).*dw_fo./sum(mob_face,2));
      zo = dt./pv(ic).*(-sgn.*(f + g.*mob(io,2)).*do_fw./sum(mob_face,2) + sgn.*g.*do_fw);
      
      J = sparse([ic1; ic; ic], [ic1; iw; io], [zd; zw; zo], nc, nc);
   end
   
   % Return function handles
   F   = @compF;
   Jac = @compJacGrav;

   % Bind prm.resred * err, @(sat) F(sat, s0, dt) and prm.lstrials
   % to function call to reduce clutter in main code.
   linsrch = @(resSol, s0, ds, dt, err, iw, io)        ...
               linesearch(resSol, ds, prm.resred * err,  ...
               @(resSol) F(resSol, s0, dt, iw, io), prm.lstrials);
end

%--------------------------------------------------------------------------

function [g, f] = getFlux(G, rock, fluid, resSol, grav_only)
% return total velocity and gravity face flux.
   is_int = all(G.faces.neighbors~=0, 2);

   [nc, nf] = deal(G.cells.num, sum(double(is_int)));
   cellNo   = rldecode(1 : nc, double(G.cells.numFaces), 2) .';

   % Indices to (internal) half-faces.
   cIntFInx = is_int(G.cellFaces(:,1));

   % Global-to-local face map for internal faces.
   %G2L      = cumsum(double(is_int));

   % Indices of internal face corresponding to cIntFInx.
   globfix  = G.cellFaces(cIntFInx, 1);

   % Renumbering globfix to 1:numel(globfix)
   %locfix   = G2L(globfix);

   % - Gravity matrix -
   % Spread gravity to cellfaces: cellFInx * g_const .* cellFace_normal

   % rho(1) - rho(2) is always correct (even when rho(2) > rho(1)), because
   % during transport (using, e.g., the 'twophaseUpwBEGrav' transport
   % solver) we only modify the *first* saturation component.
   g   = gravity() * (fluid.rho(1) - fluid.rho(2));
   dim = size(G.nodes.coords, 2);

   harm          = zeros([G.faces.num, 1]);
   renum         = zeros([G.faces.num, 1]);
   renum(is_int) = 1 : nf;

   % nKg == n' * K * g for all cellfaces.
   [K, r, c] = permTensor(rock, dim);
   nKg = sum(G.faces.normals(G.cellFaces(:,1), r) .* ...
             bsxfun(@times, K(cellNo,:), g(c)), 2);

   % Compute harmonic average of nKg on all *internal* faces.
   harm(is_int) = 1 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));

   f = resSol.faceFlux(is_int);
   f = [f;f];
   
   g = harm(is_int);
   g = [g;g];
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
   sn = resSol;
   % Geometric line search: seems pretty robust
   while fail && (i < ni),
      sn.s(:,1) = capSat(resSol.s(:,1) + pow2(ds, alph));
      sn.s(:,end) = 1-sum(sn.s(:,1:end-1), 2);
      res  = F(sn);

      alph = alph - 1;
      i    = i + 1;
      fail = ~(norm(res, inf) < target);
   end

   alph = pow2(alph + 1);      % Undo last (unneeded) scaling.
   resSol.s = sn.s;
end

%--------------------------------------------------------------------------

function nc = printif(flag, varargin)
   if flag,
      nc = fprintf(varargin{:});
   else
      nc = 0;
   end
end
