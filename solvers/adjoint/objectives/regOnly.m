function [obj] = regOnly(param, G, S, W, rock, fluid, simRes, schedule, controls, varargin)
% Regularization term only aiming to test derivative !
%
% SYNOPSIS:
%   obj = (G, S, W, rock, fluid, simRes, schedule, controls, varargin)
%
% DESCRIPTION:
%   Computes value of objective function for given simulation, and partial
%   derivatives of variables if varargin > 6
% PARAMETERS:
%   simRes      -
%
% RETURNS:
%   obj         - structure with fields
%        val    - value of objective function
%        
%   
%
%
% SEE ALSO:
%  
opt     = struct('OilPrice',                1, ...
                 'WaterProductionCost',     0.1,  ...
                 'WaterInjectionCost',      0.15,  ...
                 'RelativeDiscountFactor',  0);
opt     = merge_options(opt, varargin{:});
%-----------------------------------------------

computePartials  = (nargin > 6);
numSteps = numel(simRes);
val      = 0;
partials = repmat( struct('v', [], 'p', [], 'pi', [], 's', [], 'u', []), [numSteps 1] );

% load data from reference model, simRes_ref
load Kmodel;

m = reshape(param.m, G.cells.num, numel(schedule));
 
for step = 2 : numSteps
    
    % regularization term
    curPerm  = ( param.K ./ m );
    K        = K(:);
    regTerm  = K - curPerm ;
    numReg   = size(regTerm,1);
    covM     = eye(numReg,numReg);
    
    % Objective value: + add regularization term, free noise
    regularTerm  = regTerm'*covM*regTerm;
    
    val     = val + regularTerm;
   
    if computePartials        
        numC  = G.cells.num;
        numCF = size(G.cellFaces, 1);
        numF  = G.faces.num;
        numW  = numel(W);
        
        qw_d               = zeros(1,numW);
        partials(step).q_w = qw_d;
        
        partials(step).u   = zeros(1, numC);
        partials(step).u   = (2*regTerm ./ (m.^2))';
        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);
        
        ds                 = zeros(1, numC);
        partials(step).s   = ds;
        
        % Second order derivatives:
        d2s = zeros(numC, 1);
        partials(step).s2          = spdiags(d2s, 0, numC, numC);
        partials(step).qs          = zeros(numC, numC); % 4 is an hardcode of production well!
        partials(step).u2          = ((2*param.K ./ m.^4) + (4*regTerm ./ (m.^3)))'; 

    end
end

obj.val = val;
if computePartials, obj.partials = partials; end


function [f_w] = fracFlow(s, fluid)
% Derivative of fractional flow function of the form
%
%            s²/µw              s²                µw
%    f(s) = ---------------- = ------------ ,  mr=---
%            s²/µw+(1-s)²/µo    s²+mr*(1-s)²      µo
%    

mr   = fluid.muw/fluid.muo; % !!!!!!!!!!!!!!!!!
f_w  = ( s.^2 ) ./ (s.^2 + mr*(1-s).^2);
return


function [df] = DFracFlow(s, fluid)
% Derivative of fractional flow function of the form
%
%            s²/µw              s²                µw
%    f(s) = ---------------- = ------------ ,  mr=---
%            s²/µw+(1-s)²/µo    s²+mr*(1-s)²      µo
%
%
%            2 mr*s(1-s)
%   df(s) = ---------------
%           (s²+mr*(1-s)²)²        

mr  = fluid.muw/fluid.muo; % !!!!!!!!!!!!!!!!!
df  = ( (2*mr) * s .* (1-s) ) ./ ( (s.^2 + mr*(1-s).^2).^2 );
return

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

function [wellRates, rateSigns] = getRates(W, wellSol)
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Dw          = blkdiag( wellSys.D );
wellSigns   = ones( numel(W), 1 );
totRates    = Dw'*wellRates;
wellSigns( totRates < 0 ) = -1;
for k = 1:numel(W)
    if ~isempty(W(k).sign) % override if sign is given expl
        wellSigns(k) = W(k).sign;
    end
end

rateSigns = Dw*wellSigns;
