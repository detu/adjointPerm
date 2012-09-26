% Routines for solving transport/saturation equation.
%
% Files
%   blackoilJacobian   - Residual and Jacobian of system of nonlinear equations associated with single-point upwind scheme.
%   blackoilUpwFE      - Explicit single point upwind transport solver for Black Oil flow.
%   explicitTransport  - Explicit single point upwind transpot solver for two-phase flow.
%   findFaceMobMat     - Initialize upwind saturation index matrices for face mobility.
%   implicitBlackOil   - Implicit single point upwind transport solver for Black Oil flow.
%   implicitTransport  - Implicit single point upwind transport solver for two-phase flow.
%   inflow             - Determine upwind fluxes (and sources) from flow calculation results.
%   inflow_bo          - Determine upwind Black-Oil fluxes (and sources) from flow calucations.
%   initFaceMob        - Initialize upwind face mobility saturation indices according to Darcy flux.
%   initTransport      - Compute input to transport solver from flow calculation results.
%   newtonRaphson      - Solve non-linear equation F(z)=0 using Newton-Raphson method.
%   twophaseUpwBE      - Implicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwBEGrav  - Implicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwFE      - Explicit single point upwind solver for two-phase flow, no gravity.
%   twophaseUpwFEGrav  - Explicit single point upwind solver for two-phase flow, including gravity.
%   twophaseUpwReorder - Single point upwind solver for Buckley-Leverett flow based on reordering.

%{
#COPYRIGHT#
%}
