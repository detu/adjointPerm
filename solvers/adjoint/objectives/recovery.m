function [obj] = recovery(G, S, W, rock, fluid, simRes, schedule, controls, varargin)
% recovery -- objective function calculating water volume at last time step
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

numSteps = numel(simRes);
porv     = G.cells.volumes.*rock.poro;
obj.val  = sum( porv.*simRes(numSteps).resSol.s );

if nargin > 6
    partials = repmat( struct('v', [], 'p', [], 'pi', [], 'q_w', [], 's', [], 'u', []), [numSteps 1] );
    numCF    = size(G.cellFaces, 1);
    numC     = G.cells.num;
    numF     = G.faces.num;
    numU     = numel(controls.well);

    for k = 1 : numSteps
        ts = k;
        partials(ts).v    = zeros(1, numCF);
        partials(ts).p    = zeros(1, numC);
        partials(ts).pi   = zeros(1, numF);
        partials(ts).q_w  = zeros(1, length( vertcat(W.cells) ));
        if ts == numSteps
            partials(ts).s = porv';
        else
            partials(ts).s = zeros(1, numC);
        end
        partials(ts).u  = zeros(1, numU);
    end
    obj.partials = partials;
end

