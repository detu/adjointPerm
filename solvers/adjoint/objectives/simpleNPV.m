function [obj] = simpleNPV(G, S, W, rock, fluid, simRes, schedule, controls, varargin)
% simpleNPV - simple net-present-value function - no discount factor
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
ro      = opt.OilPrice/0.1590;
rw      = opt.WaterProductionCost/0.1590;
ri      = opt.WaterInjectionCost/0.1590;
d       = opt.RelativeDiscountFactor;
%-----------------------------------------------

computePartials  = (nargin > 6);
numSteps = numel(simRes);
val      = 0;
partials = repmat( struct('v', [], 'p', [], 'pi', [], 's', [], 'u', []), [numSteps 1] );
totTime  = max( [simRes.timeInterval] );
        
for step = 2 : numSteps
    resSol  = simRes(step).resSol;
    wellSol = simRes(step).wellSol;
    int     = simRes(step).timeInterval;
    dt      = int(2) - int(1);
    dFac    = (1+d)^(-int(2)/totTime);
    
    [wellRates, rateSigns] = getRates(W, wellSol);
    wellCells = vertcat( W.cells );
    wellSats  = resSol.s( wellCells ); 
    %f_o       = 1 - fracFlow(wellSats, fluid);  % fractional flow oil
    f_w_all   = fluid.fw(resSol);
    f_w       = f_w_all(wellCells);
    f_o       = 1 - f_w;
    injInx    = (rateSigns > 0);
    prodInx   = (rateSigns < 0);
   
    % Objective value:
   
    val   = val + dt*dFac*( - sum(  wellRates(injInx)                )*ri ...
                            - sum( -wellRates(prodInx).*f_w(prodInx) )*rw ...
                            + sum( -wellRates(prodInx).*f_o(prodInx) )*ro );
    if computePartials
        numCF    = size(G.cellFaces, 1);
        numC     = G.cells.num;
        numF     = G.faces.num;
        numU     = numel(controls.well);
        
        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);
        
        partials(step).q_w =  - dt*dFac*ri*injInx' ...
                              + dt*dFac*rw*( prodInx.*f_w )' ...
                              - dt*dFac*ro*( prodInx.*f_o )';
        Dfw_all = fluid.dfw(resSol);
        Df_w    = Dfw_all(wellCells);
        Df_o     = - Df_w;
        %Df_o  = - DFracFlow(wellSats, fluid);
        ds    = zeros(1, numC);
        ds( wellCells(prodInx) )  =  dt*dFac*rw*( wellRates(prodInx).*Df_w(prodInx) )' ...
                                    -dt*dFac*ro*( wellRates(prodInx).*Df_o(prodInx) )'; 
        partials(step).s   = ds;
        
        partials(step).u  = zeros(1, numU);
        
        % Second order derivatives:
        if isfield(fluid, 'D2fw')
            D2f_w = fluid.D2fw( wellSats );
            D2f_o = -D2f_w;
            
            d2s = zeros(numC, 1);
            d2s( wellCells(prodInx) )  =  dt*dFac*rw*( wellRates(prodInx).*D2f_w(prodInx) )' ...
                -dt*dFac*ro*( wellRates(prodInx).*D2f_o(prodInx) )';
            partials(step).s2   = spdiags(d2s, 0, numC, numC);
            % --------------------
            %IX = sparse( wellCells(prodInx), (1:nnz(prodInx))', 1, G.cells.num, nnz(prodInx));
            % ds' = IX * (dt*dFac*rw*prodRates*(IX' * Df_w(s)) + ...
            %dssl = dt*dFac*rw*( wellRates(prodInx).*D2f_w(prodInx) ) ...
            %      -dt*dFac*ro*( wellRates(prodInx).*D2f_o(prodInx) );
            %partss = IX*spdiags(dssl, 0,  nnz(prodInx), nnz(prodInx))*IX';
            
            partials(step).qs   = dt * dFac * sparse(wellCells(prodInx), find(prodInx), rw*Df_w(prodInx)-ro*Df_o(prodInx), ...
                numC, length(prodInx));
        end
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
Dw    = blkdiag( wellSys.D );
wellSigns    = ones( numel(W), 1 );
totRates = Dw'*wellRates;
wellSigns( totRates < 0 ) = -1;
for k = 1:numel(W)
    if ~isempty(W(k).sign) % override if sign is given expl
        wellSigns(k) = W(k).sign;
    end
end

rateSigns = Dw*wellSigns;
