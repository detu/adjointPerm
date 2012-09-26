function resSol = explicitBlackOil(G, resSol, wellSol, tf, pv, fluid, ...
                                   varargin)
%Explicit single point upwind transport solver for Black Oil flow.
%
% SYNOPSIS:
%    xr = explicitBlackOil(G, xr, xw, t, pv, fluid)
%    xr = explicitBlackOil(G, xr, xw, t, pv, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function blackoilUpwFE computes the solution to the Black Oil transport
%   equation
%
%        z_t + f(s)_x/B = q, with z = s / B,
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,t] while honouring a CFL time step restriction.
%
% REQUIRED PARAMETERS:
%   G      - Grid structure as described in grid_structure.
%
%   xr, xw - Reservoir and well solution structures either properly
%            initialized or the results from a previous call to flow-solver
%            function 'solveBlackOilWellSystem'.
%
%   t      - End point of time integration interval (i.e., final time)
%
%   pv     - Reservoir pore volumes.  One scalar value for each cell in the
%            model.
%
%   fluid  - Black Oil fluid object as defined by, e.g., function
%            'initBlackoilFluid'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - Well structure as defined by function 'addWell'.  May be
%            empty (i.e., wells = [], default value) which is interpreted
%            as a model without any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = [], default
%            value) which is interpreted as all external no-flow
%            (homogeneous Neumann) conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = [], default value)
%            which is interpreted as a reservoir model without explicit
%            sources.
%
% RETURNS:
%   xr - Updated reservoir solution structure with new values for the
%        fields:
%          - s -- Array, size G.cells.num-by-np, of phase saturations.  One
%                 column for each fluid phase.
%          - z -- Array, size G.cells.num-by-np, of phase surface volumes.
%                 One column for each fluid phase.
%
% NOTE:
%   This function does not presently honour a CFL-type condition on the
%   time step size.
%
% EXAMPLE:
%   blackOilTube
%
% SEE ALSO:
%   initBlackoilFluid, blackoilJacobian, implicitBlackOil.

%{
#COPYRIGHT#
%}

% $Date: 2009-06-05 19:19:30 +0200 (fr, 05 jun 2009) $
% $Revision: 2338 $

   opt = struct('verbose', true, ... % Emit progress reports
                'bc',      [],   ...
                'src',     [],   ...
                'wells',   []);

   opt = merge_options(opt, varargin{:});

   F = blackoilJacobian2(G, resSol, wellSol, pv, fluid, ...
                         opt.wells, opt.src, opt.bc);

                    
   cellNo = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
   cFlux  = faceFlux2cellFlux(G, resSol.faceFlux);
   fluxp  = accumarray(cellNo, max(cFlux,  0));
   fluxn  = accumarray(cellNo,-min(cFlux,  0));
   tau    = min(pv./max([fluxp, fluxn], [], 2));
   
   z  = reshape(resSol.z', [], 1);
   dt = min(0.25*tau, tf);
   t  = 0;

   while t < tf,
      z(:) = z - F(z, z, dt);

      t    = t + dt;
      dt   = min(dt, tf - t);
   end

   resSol.z  = reshape(z, size(resSol.z'))';
   [u,u,u,u] = fluid.pvt(resSol.cellPressure, resSol.z);
   resSol.s  = bsxfun(@rdivide, u, sum(u,2));
end
