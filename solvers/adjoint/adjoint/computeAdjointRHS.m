function b = computeAdjointRHS(G, W, f_res) 
% Compute adjoint 'pressure' rhs 
% for use in e.g. function solveAdjointPressureSystem
%
% SYNOPSIS:
%   b = computeAdjoint(G, W, f_res, f_w)
% 
% DESCRIPTION:
% Computes adjoint 'pressure' rhs as input to solveIncompFlow 
% in function solveAdjointPressureSystem
%
% PARAMETERS:
%   G     - Grid data structure.
%  
%   W     - Well structure as defined by addWell &c.
%
%   f_res - Adjoint reservoir 'pressure' condtions
%
% RETURNS:
%   b     - Ajoint pressure rhs to be passed directly as option
%           'rhs' to solveIncompFlow.
%

assert( numel(f_res) == size(G.cellFaces,1)); 
f_w = [];

% unpack adjoint well 'pressure' conditions
if ~isempty(W),  
   S   = [ W.S   ];
   RHS = [ S.RHS ];
   f_w = {RHS.f};
   f_w = vertcat(f_w{:});
   % start of changes by EKA
   h_w = zeros(numel(vertcat(W.val)),1);  % BHP-controlled
%    h_w = {RHS.h};
%    h_w = vertcat(h_w{:});
   % end of changes by EKA
end

% b = [f_res; f_w; g; h_res; h_w]
b    = cell([1, 3]);
b{1} = vertcat(f_res, f_w);
b{2} =         zeros([G.cells.num, 1]);
b{3} = vertcat(zeros([G.faces.num, 1]), h_w);
end
