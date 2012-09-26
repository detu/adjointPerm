function [obj] = PermMismatch(param, G, S, W, rock, fluid, simRes, schedule, controls, varargin)
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
%-----------------------------------------------

computePartials  = (nargin > 7);
numSteps = numel(simRes);
val      = 0;
partials = repmat( struct('v', [], 'p', [], 'pi', [], 's', [], 'u', []), [numSteps 1] );

% load data from reference model, simRes_ref
load simResSmallRef;
load Kmodel;

m = param.m;

for step = 2 : numSteps
    
    % regularization term added only at final time
    if step == numSteps
        % regularization term
        curPerm  = ( param.K ./ m );
        K        = K(:);
        regTerm  = K - curPerm ;
        numReg   = size(regTerm,1);
        covM     = eye(numReg,numReg);
        regularTerm  = regTerm'*covM*regTerm;
        
        val = val + regularTerm;
    end
    
    if computePartials        
        numC  = G.cells.num;
        numCF = size(G.cellFaces, 1);
        numF  = G.faces.num;
        numW  = numel(W);
        
        qw_d               = zeros(1,numW);
        partials(step).q_w = qw_d;
        
        partials(step).v   = zeros(1, numCF);
        dp                 = zeros(1, numC);
        partials(step).p   = dp;
        partials(step).pi  = zeros(1, numF);
        ds                 = zeros(1, numC);
        partials(step).s   = ds;
        
        if step == numSteps
            partials(step).u   = (2*regTerm ./ (m.^2))';
        else
            partials(step).u   = zeros(1, numC);
        end
        
    end
end

obj.val = val;
if computePartials, obj.partials = partials; end

