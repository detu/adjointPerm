%
%  function [ Hv ] = Hessvec( v, u, usr_par) 
%
%     Purpose:
%
%     Compute the product of the Hessian objective function f(y(u),u)
%     times a vector v, where f(y(u),u) is the objective function 
%     in the optimal control of the unsteady Burger's equation. 
%
%
%     Parameters
%
%     On entry:
%
%     v      Direction v.
%            v((i-1)*(nx+1)+1), ..., v(i*(nx+1))
%            direction at time (i-1)*Deltat, i = 1, ..., nt+1
%
%     u      Control  u.
%            u((i-1)*(nx+1)+1), ..., u(i*(nx+1))
%            controls at time (i-1)*Deltat, i = 1, ..., nt+1
%
%     usr_par user defined parameter. Used to pass problem
%            specific information.
% 
%     On return:
%
%     Hv      Value of Hessian times vector product.
%
%
% Version June 6, 2008
% Modified by Eka May 7, 2009
% Changing for Oil Reservoir Model
% Matthias Heinkenschloss
%
function [ Hv ] = hess2( dv, u, usr_par) 

[simRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction] = deal(usr_par{:});

% forward run to solve derivative of the states / w
% pressure equation and its derivatives : w_p, w_v, w_phi

% Objective function NPV
controls  = updateControls(controls, dv);
schedule  = updateSchedule(controls, schedule);

% FORWARD SOLVE
simRes    = runSchedule(simRes(1).resSol, G, S, W, rock, fluid, schedule, 'VerboseLevel', 0);
numSteps  = numel(schedule);

% ADJOINT SOLVE
adjRes    = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls,objectiveFunction);

dv        = reshape(dv,numSteps,numSteps);
simResDer = [];

for curStep = 1 : numSteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE EQUATION
    dt      = simRes(curStep+1).timeInterval * [-1 1]';
    
    % Generate right-hand-side
    numCF   = size(S.B, 1);
    numW    = numel(W);
    dim     = fluid.DLtInv(simRes(curStep).resSol.sw);
    numC    = G.cells.num;
    
    % Update B^{n-1}v^{n-1} part
    v       = simRes(curStep).resSol.cellFlux;   % v^{n-1}
    if curStep == 1
        w_s = zeros( numC, 1 );
    else
        w_s = simResDer(curStep).resSol.sw;   % w_s^{n-1}
    end
    RHS     = -S.B*spdiags( (S.C*dim).*v , 0, numCF, numCF)* S.C * w_s;
    
    % Update B_w^{n-1}q_w^{n-1} part
    for wellNr = 1:numW
        w   = W(wellNr);
        if curStep == 1
            q   = zeros( size(simRes(curStep+1).wellSol(wellNr).flux) ); % q_w^{n-1}
        else
            q   = simRes(curStep).wellSol(wellNr).flux; % q_w^{n-1}
        end
        W(wellNr).S.RHS.f =  w.S.B * diag(( ( (w.S.C*dim).*q ) )) * w.S.C * w_s;
    end
    
    % should be checked, B_w*q_w goes to reservoir or well equation ?!
    S.RHS.f_bc = RHS;  
    
    % Set f and h in W(i).S.RHS equal u' and zero. Leave naumannFaces and dirichletFaces as is 
%     inx = 0;
    for wellNr = 1 : numel(W)
%         numCells = length( W(wellNr).cells );
%         W(wellNr).S.RHS.f = obj.partials(curStep).q_w( inx+1 : inx+numCells )';
        W(wellNr).S.RHS.f = W(wellNr).S.RHS.f + dv(wellNr,curStep);
        W(wellNr).S.RHS.h = zeros( size(W(wellNr).S.RHS.h) );
%         inx = inx + numCells;
    end

    % Solve linear system based on s^{n-1}
    [resSol, wellSol] = solveWellSystem(simRes(curStep).resSol, G, S, W, fluid);

    % Update adjRes !!! Note minuses in front of pressure and wellrates in
    % forward system, but not in adjoint, thus set minus here
    simResDer(curStep+1).resSol.cellFlux     = resSol.cellFlux;
    simResDer(curStep+1).resSol.cellPressure = - resSol.cellPressure;           % !!!minus
    simResDer(curStep+1).resSol.facePressure = resSol.facePressure;
    simResDer(curStep+1).wellSol             = wellSol;
    for k = 1 : numel(simResDer(curStep+1).wellSol)
        simResDer(curStep+1).wellSol(k).flux = - simResDer(curStep+1).wellSol(k).flux;     % !!!minus
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SATURATION EQUATION
    
    % Generate system matrix
    numC    = G.cells.num;
    PV      = G.cells.volumes.*rock.poros;
    invDPV  = spdiags(1./PV, 0, numC, numC);
    DDf     = spdiags( fluid.Dfw(simRes(curStep+1).resSol.sw), 0, numC, numC);
    At      = generateUpstreamTransportMatrix(G, S, W, simRes(curStep+1).resSol, ...
              simRes(curStep+1).wellSol, 'Transpose', false);
    systMat = speye(numC, numC) - dt * ( DDf * At * invDPV);   % system matrix
    clear PV invDPV DDf At
    
    % Generate RHS, that is f-part (rest is zero)
    PV      = G.cells.volumes.*rock.poros;
    invPV   = 1./PV;
    f_w     = fluid.fw( simRes(curStep+1).resSol.sw );
    if curStep == 1
        w_v = zeros( size(simResDer(curStep+1).resSol.cellFlux) );
    else
        w_v     = simResDer(curStep+1).resSol.cellFlux;        % w_v^n
    end
    
    % Flux-matrix: A.i, A.j, A.qMinus
    [A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, simRes(curStep+1).resSol, ...
                         simRes(curStep+1).wellSol, 'VectorOutput', true);
    dQPluss =  double( signQ > 0 );
    dQMinus = -double( signQ < 0 );

%     RHS    = dt*( f_w(A.j).*( S.C*(invPV.*w_v) ) ) ...
%             + dt*( S.C*( (-f_w.*dQMinus + dQPluss).*( invPV.*w_v ) ) );
    RHS    = dt* ( diag( f_w(A.j) )* S.C * diag(invPV) )' * w_v ...
            + dt* ( S.C* diag( (-f_w.*dQMinus + dQPluss).*( invPV ) ))'* w_v   ;
    
    if curStep > 1
        RHS = RHS + simResDer(curStep).resSol.sw;
    end
    simResDer(curStep+1).resSol.sw  = systMat \ RHS;
    
end

adjResDer= [];
obj      = objectiveFunction(G, S, W, rock, fluid, simRes, schedule, controls);
for k = numSteps : -1 : 1

    W    = updateWells(W, schedule(k));
    
     % ---- Analogue of Saturation Equation 
    adjResDer = solveDATransportSystem(G, S, W, rock, fluid, simRes, adjResDer, obj, simResDer);
    
    % ---- Analogue of Pressure Equation
    adjResDer = solveDAPressureSystem(G, S, W, rock, fluid, simRes, adjResDer, obj);

end

% Finally, compute the Hessian !
hess = computeHessian(W, adjResDer, schedule, controls);
Hv    = cell2mat( {hess{:}}');



% End of hessvec.

function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end
