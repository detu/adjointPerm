function [obj] = waterCutEnKF(G, W, fluid, g_mEns, g_m, asimStep, cov, pPerm, param, simRes, varargin)
% TO DO: - GANTI TOTAL FUNCTION INI, gak perlu iterate setiap numStep,
% karena pake local KF, numstep jadi input ke function ini
%        - REFENCE WATER CUT input juga dari g_mEns
%
%
% Water Cut Match Objective Function + Regularization term
%
% SYNOPSIS:
%   obj = (G, S, W, rock, fluid, simRes, schedule, controls, varargin)
%
% DESCRIPTION:
%   Computes value of objective function for given simulation, and partial
%   derivatives of variables if varargin > 7
% PARAMETERS:
%   simRes      -
%
% RETURNS:
%   obj         - structure with fields
%        val    - value of objective function
%        
%
% SEE ALSO:
%  

computePartials  = (nargin > 9);
numSteps = 2;
partials = repmat( struct('v', [], 'p', [], 'pi', [], 's', [], 'u', []), [numSteps 1] );

K = param.K;
m = param.m; %prior permeability

% construct objective function given input time step
% 1. mismatch term + measurement noise
measureNoise = 0.05*g_m(:,asimStep);
misMatch     = g_mEns - (g_m(:,asimStep) + measureNoise);

% 2. regularization term, i.e., current permeability with its prior
curPerm  = ( K ./ m );
priPerm  = ( K ./ pPerm );
curPerm  = curPerm(:);
priPerm  = priPerm(:);
regTerm  = priPerm - curPerm;

% 3. include covariance information
% mismatchTerm = (misMatch'/cov.Cd)*misMatch;
% regularTerm  = (regTerm'/ cov.Cm)*regTerm;
% val          = 1e-16*mismatchTerm + 1e-5*regularTerm;
Cd           = inv(cov.Cd);
Cm           = inv(cov.Cm);
mismatchTerm = misMatch'*Cd*misMatch;
regularTerm  = regTerm'*Cm*regTerm;
wm           = -1*order(mismatchTerm);
wr           = -1*order(regularTerm);
wm           = 10^(wm);
wr           = 10^(wr);
val          = wm*mismatchTerm + wr*regularTerm;

step = numSteps;

if computePartials
    resSol  = simRes(asimStep+1).resSol;
    wellSol = simRes(asimStep+1).wellSol;
    
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells );
    prodInx   = (rateSigns < 0);
    
    numC  = G.cells.num;
    numCF = size(G.cellFaces, 1);
    numF  = G.faces.num;
    numW  = numel(W);
    
    qw_d  = zeros(1,numW);
    partials(step).q_w = qw_d;
    partials(step).u   = (2*wr*Cm*regTerm ./ (m.^2))';
    %partials(step).u   = (2*wr*regTerm ./ (m.^2))';
    
    partials(step).v   = zeros(1, numCF);
    partials(step).p   = zeros(1, numC);
    partials(step).pi  = zeros(1, numF);
    
    
    Dfw_all                    = fluid.dfw(resSol);
    Df_w                       = Dfw_all(wellCells(prodInx));
    ds                         = zeros(1, numC);
    ds( wellCells(prodInx) )   = 2*wm*Cd*misMatch .* Df_w;
    %ds( wellCells(prodInx) )   = 2*wm*misMatch .* Df_w;
    partials(step).s           = ds;
    
    % Second order derivatives:
    D2f_w                      = fluid.d2fw(struct('s', wellSats) );
    d2s                        = zeros(numC, 1);
    d2s( wellCells(prodInx) )  = 2*wm*Cd* ( -Df_w.*Df_w +  misMatch .* D2f_w(prodInx) );
    partials(step).s2          = spdiags(d2s, 0, numC, numC);
    partials(step).qs          = zeros(numC, length(prodInx));
    partials(step).u2          = 2*wr*Cm* ( ( param.K ./ m.^4 ) +  ( -2*regTerm ./ (m.^3 ) ) );
    
end

obj.val = val;
if computePartials, obj.partials = partials; end

return;
