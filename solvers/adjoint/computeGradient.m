function grad = computeGradient(W, adjRes, schedule, controls)
% compute gradient for control variables and project according to
% linEqConst A*u = b  => A*grad = 0
% Thus the projected gradient is 
%   grad_p  = (I - A'*inv(A*A')*A)*grad 


bhpWells  = find( strcmp('bhp', {W.type}) );
rateWells = find( strcmp('rate', {W.type}) );
S         = [ W.S ];
Dw        = blkdiag( S.D );
DwD       = Dw(:, bhpWells);

% collect all linEqConst
ec = controls.linEqConst;
if ~isempty(ec)
    A = ec.A;
    projector = eye( size(A, 2) ) - A'*( (A*A')\A );
end

[A_N, b_N, A_D, b_D] = controls2Wells(W, schedule, controls);

for k = 1 : numel(A_N)
    adjWellPres = [adjRes(k+1).wellSol.pressure]';
    l_p         = adjWellPres(rateWells);
    l_q         = vertcat(adjRes(k+1).wellSol.flux);
        
    grad_k  =  A_N{k}'*l_p + A_D{k}'*DwD'*l_q;   % non-projected gradient
    if ~isempty(ec) 
        grad{k}  = projector * grad_k;
    else
        grad{k}  = grad_k;
    end
end

if controls.numControlSteps == 1
    gradMat = cell2mat(grad);
    avgVal  = mean(gradMat, 2);
    for k = 1 : numel(A_N)
        grad{k} = avgVal;
    end
end

    