function S = newtonRaphson(S, dT, sources, gm, porevolumes, varargin)
% newtonRaphson -- Solves the water saturation equation with
%                  Newton-Raphson method with a fixed sub time step of
%                  dT/2^k.
%
% SYNOPSIS:
%   S = newtonRaphson(S, dT, q, gm, porvol)
%   S = newtonRaphson(S, dT, q, gm, porvol, 'pn', pv, ...)
%
% REQUIRED PARAMETERS:
%   s        - Initial water saturation (at time=0).  One scalar value for
%              each cell in the grid.
%
%   t        - Length of time step.
%
%   q        - Accumulated sources (typically contributions from wells).  One
%              scalar value for each cell in the grid.
%
%   gm       - Upwind/inflow matrix of fluxes into each cell, created (e.g.)
%              by function INFLOW. gm(i,j) is the (positive) flux from cell
%              j to cell i.
%
%   porvol   - Reservoir pore volumes, one scalar value for each cell in the
%              grid.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   verbose  - Whether or not time integration progress should be reported to
%              the screen.
%              Default value: verbose = false.
%
%   mr       - Mobility ratio, i.e., ratio of water mobility to oil mobility.
%              Default value: mr = 1.0.
%
%   nltol    - Absolute tolerance of iteration.  The numerical solution must
%              satisfy the condition
%
%                NORM(S-S0 + dt/porvol(out - in) - Q, INF) <= nltol
%
%              at all times in the interval [0,t].
%
%              Default value: nltol = 1.0e-3.
%
%   maxnewt  - Maximum number of inner iterations in Newton-Raphson method.
%              Default value: maxnewt = 25.
%
% RETURNS:
%   s        - Water saturation, one scalar value for each cell, at time t.
%
% NOTE:
%   Original code by J.E. Aarnes.
%
% SEE ALSO:
%   INFLOW.

% $Id: newtonRaphsonJEA.m 1264 2009-01-21 13:55:46Z jrn $

prm  = struct('verbose',  false,  ...  % emit progress reports
              'mr',       1.0,    ...  % mobility ratio
              'nltol',    1.0e-3, ...  % non-linear residual tolerance
              'maxnewt',  25);         % max no. of NR iterations

prm = merge_options(prm, varargin{:});

% Support multi-D shaped saturation input array.
sz  = size(S);
S   = reshape(S, [], 1);

N   = numel(S);
dv  = prm.mr;  % µw/µo
gp  = spdiags(-sum(gm, 1)', 0, N, N);
F   = -(gp + gm);

dsnorm = 1;

IT  = 0;
S00 = S;

while dsnorm > prm.nltol, %~convergence
  PV = pow2(dT, -IT) ./ porevolumes;
  B  = spdiags(PV, 0, N, N) * F;
  fi = PV .* sources; I=0;

  % Try using 2^IT substeps
  while I < 2 ^ IT,

    S0     = S;
    dsnorm = 1;
    it     = 0;
    I      = I + 1;

    while (dsnorm > prm.nltol) && (it < prm.maxnewt),
      % evaluate flux function (fw), and its derivative fd
      fn = S.^2 + dv*(1-S).^2; ft = (2*dv) * S .* (1 - S);
      fw = S.^2 ./ fn;         fd = ft ./ fn.^2; %fw

      % Jacobian of G
      GG = speye(N) - B*spdiags(fd, 0, N, N);

      % G
      gg = S0-S + fi + B*fw;
      ds = GG \ gg;

      % Add (capped) increment ds (to yield 0 <= S <= 1 always)
      S1  = min(S + ds, 1); ds = S1 - S;
      S   = S + ds;
      dsnorm = norm(ds);
      it = it + 1;
    end

    if dsnorm > prm.nltol
      % Set I to 2^IT to finish 'while I < 2^IT...' loop
      I = 2 ^ IT;
      % Set S to S at start of time step (i.e. redo time step)
      S = S00;
    end
  end

  if dsnorm > prm.nltol,
    % Cut time step in half
    IT = IT + 1;
  end
end

S = reshape(max(S, 0), sz);
