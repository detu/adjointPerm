function [obj] = Jac_Swn(G, S, W, rock, fluid, simRes, schedule, controls, varargin)
% Jacobian of water saturation (as state constraint) constrained case
% MUST BE TESTED BEFORE USE THIS BY RUNNING COMPAREGRADIENT !!!
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
% opt     = struct('OilPrice',                1, ...
%                  'WaterProductionCost',     .15,  ...
%                  'WaterInjectionCost',      .1,  ...
%                  'RelativeDiscountFactor',  0);
opt     = struct('OilPrice',                1, ...
                 'WaterProductionCost',     0,  ...
                 'WaterInjectionCost',      0,  ...
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


% water saturation bound
Sbound = 0.65;

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
    if step == numSteps
        
        % constrained water saturation at final time
        val = sum( max(( wellSats(prodInx) - Sbound ), 0)  );
        
    end
    

    if computePartials
        numCF    = size(G.cellFaces, 1);
        numC     = G.cells.num;
        numF     = G.faces.num;
        numU     = numel(controls.well);
        
        partials(step).v   = zeros(1, numCF);
        partials(step).p   = zeros(1, numC);
        partials(step).pi  = zeros(1, numF);
        
        partials(step).q_w =  zeros(size(wellCells))';
        
        ds    = zeros(1, numC);

        if step == numSteps
            sw = wellSats(prodInx);
            wc = wellCells(prodInx);
            for ind = 1: numel( wellCells(prodInx) )
                if ( sw(ind) ) > Sbound
                    ds( wc(ind) ) = 1;
                else
                    ds( wc(ind) ) = 0;
                end
            end
                
%             if wellSats(prodInx) > Sbound
%                 ds( wellCells(prodInx) )  = 1;
%             else
%                 ds( wellCells(prodInx) )  = 0;
%             end
        end
        partials(step).s   = ds;
        
        partials(step).u  = zeros(1, numU);
    end
end

obj.val = val;
if computePartials, obj.partials = partials; end


function [f_w] = fracFlow(s, fluid)
% Derivative of fractional flow function of the form
%
%            s??/??w              s??                ??w
%    f(s) = ---------------- = ------------ ,  mr=---
%            s??/??w+(1-s)??/??o    s??+mr*(1-s)??      ??o
%    

mr   = fluid.muw/fluid.muo; % !!!!!!!!!!!!!!!!!
f_w  = ( s.^2 ) ./ (s.^2 + mr*(1-s).^2);
return


function [df] = DFracFlow(s, fluid)
% Derivative of fractional flow function of the form
%
%            s??/??w              s??                ??w
%    f(s) = ---------------- = ------------ ,  mr=---
%            s??/??w+(1-s)??/??o    s??+mr*(1-s)??      ??o
%
%
%            2 mr*s(1-s)
%   df(s) = ---------------
%           (s??+mr*(1-s)??)??        

mr  = fluid.muw/fluid.muo; % !!!!!!!!!!!!!!!!!
df  = ( (2*mr) * s .* (1-s) ) ./ ( (s.^2 + mr*(1-s).^2).^2 );
return


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
