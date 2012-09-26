function [linSimRes] = runLinearized(simRes, G, S, W, rock, fluid, schedule, varargin)
% runSchedule -- Run simulation based on schedule.
%
% SYNOPSIS:
%   simRes = runSchedule(resSolInit, G, S, W, fluid, schedule, pn, pv, ...)
%
% DESCRIPTION:
%   
% PARAMETERS:
%   resSolInit  -
%   G           - Grid data structure.
%   S           -
%   W           -
%   fluid       -
%   schedule    -
%
%
% RETURNS:
%   simRes      - (numSteps+1) x 1 structure having fields
%                   - timeInterval
%                   - resSol
%                   - wellSol
%
%
% SEE ALSO:
%  
opt     = struct('Verbose',  false , ...
                 'VerboseLevel', 2);
opt     = merge_options(opt, varargin{:});
verboseLevel2 = opt.Verbose || (opt.VerboseLevel == 2);
verboseLevel1 = opt.Verbose || (opt.VerboseLevel > 0);

numSteps = numel(schedule);
resSol   =  initResSol(G, 0);

numCF    = size(S.B, 1);
numW     = numel(W);
%pv       = G.cells.volumes.*rock.poros;
%mobRatio = (fluid.muw/fluid.muo);  % !!!!!! Note this is inverse of usual def.

% Initial conditions
linSimRes(1).timeInterval  = [0 0];
linSimRes(1).resSol        = resSol;
linSimRes(1).wellSol       = [];
if verboseLevel2, dispSchedule(schedule); end

%dim = fluid.DLtInv(simRes(curStep).resSol.sw);

if verboseLevel1, fprintf('\n******* Starting linearized forward simualtion *******\n'); end
for k = 1 : numSteps
    if verboseLevel1, fprintf('Time step %3d of %3d,   ', k, numSteps); end
    W        = updateWells(W, schedule(k));
    interval = schedule(k).timeInterval;
    dt       = interval(2) - interval(1);
    dim      = fluid.DLtInv(simRes(k).resSol.sw);
    %f_wp1    = fluid.fw( simRes(k+1).resSol.sw );
    
    % ---- Pressure Equation -----
    % Include partial derivatives wrt s^{n-1} on RHS
    v = simRes(k+1).resSol.cellFlux;
    S.RHS.f_bc = - S.B*spdiags( (S.C*dim).*v , 0, numCF, numCF)* S.C * linSimRes(k).resSol.sw;
    
    % Update B_w^{n-1}q_w^{n-1} part
   
    for wellNr = 1:numW
        w   = W(wellNr);
        q = simRes(k+1).wellSol(wellNr).flux; 
        %W(wellNr).S.RHS.f =  W(wellNr).S.RHS.f + w.S.B * diag(( ( (w.S.C*dim).*q ) )) * w.S.C * linSimRes(k).resSol.sw;
        W(wellNr).S.RHS.f =  W(wellNr).S.RHS.f ...
            + w.S.B * spdiags(( ( (w.S.C*dim).*q ) )) * w.S.C * linSimRes(k).resSol.sw;
    end
    
    if verboseLevel1, fprintf('Pressure:'); tic; end
    [resSol, wellSol] = solveMixedWellSystem(simRes(k).resSol, G, S, W, fluid, 'Verbose', verboseLevel2);
    if verboseLevel1, t = toc; fprintf('%9.3f sec,   ', t); end
    
        
    % ---- Saturation Equation ---
    if verboseLevel1, fprintf('Transport:'); tic; end
    
    % Generate system matrix
    numC    = G.cells.num;
    PV      = G.cells.volumes.*rock.poros;
    invPV   = 1./PV;
%     invDPV  = spdiags(1./PV, 0, numC, numC);
    invDPV  = spdiags(invPV , 0, numC, numC);
    DDf     = spdiags( fluid.Dfw(simRes(k+1).resSol.sw), 0, numC, numC);
    [A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(k+1).resSol, ...
                         simRes(k+1).wellSol, 'VectorOutput', true);
                     
    AMat = sparse(A.i, A.j, -simRes(k+1).resSol.cellFlux, numC, numC) + spdiags(A.qMinus, 0, numC, numC); 
    systMat = speye(numC, numC) - dt * ( invDPV * AMat * DDf);   % system matrix
    
    dQPluss =  double( signQ > 0 );
    dQMinus = -double( signQ < 0 );
    
    f_w     = fluid.fw( simRes(k+1).resSol.sw );
    
%     RHS    = - dt* ( diag( f_w(A.j) )* S.C * diag(invPV) )' * resSol.cellFlux ...
%              + dt* ( S.C* diag( (-f_w.*dQMinus + dQPluss).*( invPV ) ))'* resSol.cellFlux ...
%              + linSimRes(k).resSol.sw;

    RHS    = - dt* ( spdiags( f_w(A.j) , 0, numCF, numCF)* S.C * invDPV )' * resSol.cellFlux ...
             + dt* ( S.C* spdiags( (-f_w.*dQMinus + dQPluss).*( invPV ), 0, numC, numC ) )'* resSol.cellFlux ...
             + linSimRes(k).resSol.sw ;
    
    resSol.sw  = systMat \ RHS;
    clear PV invDPV DDf At
    
    if verboseLevel1, t = toc; fprintf('%9.3f sec\n', t); end                         
    
    % update simRes structure
    linSimRes(k+1).timeInterval  = interval;
    linSimRes(k+1).resSol        = resSol;
    linSimRes(k+1).wellSol       = wellSol;
end

    