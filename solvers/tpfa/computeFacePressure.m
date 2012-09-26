function xr = computeFacePressure(xr, G, T, fluid, varargin)
%Compute face pressure using two-point flux approximation.
%
% SYNOPSIS:
%   xr = computeFacePressure(xr, G, T, fluid)
%   xr = computeFacePressure(xr, G, T, fluid, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   xr -     Reservoir solution structures either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompTPFA' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeTrans'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the field:
%          - facePressure -- Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%
%
% SEE ALSO:
%   computeTrans, addBC, addSource, initSingleFluid, initResSol,
%   incompTPFA.

%{
#COPYRIGHT#
%}

% $Date: 2009-10-06 10:29:14 +0200 (ti, 06 okt 2009) $
% $Revision: 2955 $

   opt = struct('bc', [], 'src', []);
   opt = merge_options(opt, varargin{:});
   
   [c, rho, mu, u] = fluid.pvt(xr.cellPressure, xr.z);
   
   mob    = fluid.relperm(bsxfun(@rdivide, u, sum(u, 2)))./mu;
   totmob = sum(mob, 2);
   omega  = bsxfun(@rdivide, rho.*mob, totmob);
   
   [ff, gg, hh, grav, dF, dC] = computePressureRHS(G, omega, ...
                                                   opt.bc, opt.src); %#ok
   % Reconstruct face pressures and fluxes.
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   
   
   % Face transmissibility = harmonic average of half-transmissibilities
   
   T      = T.*totmob(cellNo);
   %ft     = 1./accumarray(G.cellFaces(:,1), 1./T, [G.faces.num, 1]);

   
   p  = xr.cellPressure;
   fp =  ...
          accumarray(G.cellFaces(:,1), (p(cellNo)+grav).*T, [G.faces.num,1])./ ...
          accumarray(G.cellFaces(:,1), T, [G.faces.num,1]);

   %{
   % Neumann faces
   b     = any(G.faces.neighbors==0, 2);
   fp(b) = fp(b) - hh(b)./ft(b);


   % Dirichlet faces
   fp(dF) = dC;
   %}
   i = ~any(G.faces.neighbors==0, 2);
   xr.facePressure(i) = fp(i);
