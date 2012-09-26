% Routines supporting the mimetic method for the pressure equation.
%
% Files
%   computeMimeticIP            - Compute mimetic inner product matrices.
%   assembleWellSystem          - Generate pressure linear system components for wells.
%   packageWellSol              - Convert well fluxes and pressures to structure form.
%   solveBlackOilWellSystem     - One iteration of succsubst algorithm on Black Oil pressure system.
%   solveIncompFlow             - Solve incompressible flow problem (fluxes/pressures).
%   succSubst                   - Successive substitution algorithm for compressible pressure system.
%   tpf                         - Solve symmetric system of linear eqns using reduction to a mixed system.
%   unpackWellSystemComponents  - Extract hybrid linear system components from wells.

%{
#COPYRIGHT#
%}
