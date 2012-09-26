function [obj] = waterCutReg(param, G, S, W, rock, fluid, simRes, schedule, controls, varargin)
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

computePartials  = (nargin > 7);
numSteps = numel(simRes);
val      = 0;
partials = repmat( struct('v', [], 'p', [], 'pi', [], 's', [], 'u', []), [numSteps 1] );

% load data from reference model, simRes_ref, covariance matrices
load simResSmallRef;
load Kmodel;
load covMats;

m = param.m;
 
for step = 2 : numSteps
    % model
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;
    
    % measurement/reference
    resSolref  = simRes_refSmall(step).resSol;
    wellSolref = simRes_refSmall(step).wellSol;
    
    
    % model
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells ); 

    
    % measurement/reference
    [wellRates_ref, rateSigns_ref] = getRates(W, wellSolref);

    % model
    f_w_all   = fluid.fw(resSol);
    f_w       = f_w_all(wellCells);
    f_o       = 1 - f_w;
    injInx    = (rateSigns > 0);
    prodInx   = (rateSigns < 0);
    
    % measurement/reference
    f_w_all_ref   = fluid.fw(resSolref);
    f_w_ref       = f_w_all_ref(wellCells);
    f_o_ref       = 1 - f_w_ref;
    injInx_ref    = (rateSigns_ref > 0);
    prodInx_ref   = (rateSigns_ref < 0);
    
    % model - liquid rate at producer wells
    g_m     = f_w( prodInx ) ;
    
    % measurement - liquid rate at producer wells
    d_obs    = f_w_ref( prodInx_ref );
   
    misMatch = d_obs - g_m;
    numDiff  = size(misMatch,1);
    covD     = eye(numDiff,numDiff);
    mismatchTerm = misMatch'*covD*misMatch;
    
    wm  = 1;
%     wm  = 1e3;
    val = val + wm*mismatchTerm;
    
    % regularization term added only at final time
    if step == numSteps
        % regularization term
        curPerm  = ( param.K ./ m );
        K        = K(:);
        regTerm  = K - curPerm ;
        numReg   = size(regTerm,1);
        covM     = eye(numReg,numReg);
        regularTerm  = regTerm'*covM*regTerm;
        
        wr  = 1;
        val = val + wr*regularTerm;
    end
   
    if computePartials        
        numC  = G.cells.num;
        numCF = size(G.cellFaces, 1);
        numF  = G.faces.num;
        numW  = numel(W);
        
        qw_d               = zeros(1,numW);
        partials(step).q_w = qw_d;
        
        if step == numSteps
            partials(step).u   = (2*wr*regTerm ./ (m.^2))';
        else
            partials(step).u   = zeros(1, numC);
        end
        
        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);
        

        Dfw_all                    = fluid.dfw(resSol);
        Df_w                       = Dfw_all(wellCells(prodInx));
        ds                         = zeros(1, numC);
%         ds( wellCells(prodInx) )   = 2*wm*Df_w .* misMatch ;
        ds( wellCells(prodInx) )   = 2*wm*misMatch .* Df_w;
        partials(step).s           = ds;
        
        % Second order derivatives:
        D2f_w                      = fluid.d2fw(struct('s', wellSats) );
        d2s                        = zeros(numC, 1);
%         d2s( wellCells(prodInx) )  = 2*wm* misMatch .* ( D2f_w(prodInx) + Df_w ); % CHECK AGAIN !!!
        d2s( wellCells(prodInx) )  = 2*wm* ( -Df_w.*Df_w +  misMatch .* D2f_w(prodInx) ); 
        partials(step).s2          = spdiags(d2s, 0, numC, numC);
        partials(step).qs          = zeros(numC, length(prodInx));
        
        if step == numSteps
            partials(step).u2      =  2*wr* ( ( param.K ./ m.^4 ) +  ( -2*regTerm ./ (m.^3 ) ) );
        end

    end
end

obj.val = val;
if computePartials, obj.partials = partials; end


%{
function [wellRates, rateSigns] = getRates(W, wellSol)
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Dw    = blkdiag( wellSys.D );
if isfield(W, 'sign')
    wellSigns = vertcat(W.sign);
else
    wellSigns    = ones( numel(W), 1 );
    totRates = Dw'*wellRates;
    wellSigns( totRates < 0 ) = -1;
end
%}


