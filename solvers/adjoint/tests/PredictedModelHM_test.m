% PredictedModel HM
% Initial Predicted Model to match SimpleModelHM
% a 4 x 4 reservoir for history matching

nx = 4; ny = 4; nz = 1;
G  = cartGrid([nx ny nz]);
G  = computeGeometry(G);

% define K -> K/m, where m is permeability multiplier
% let K = one
m     = [ 2 2 3 3; ...
          2 2 3 3; ...
          4 4 1 1; ...
          4 4 1 1 ];
Kreal = [ 1 1 1 1; ...
          1 1 1 1; ...
          1 1 1 1; ...
          1 1 1 1 ];
K     = Kreal ./ m;

perm  = reshape(K, [1,G.cells.num]);
rock.perm = perm'*100*milli*darcy;

% load permSmallInit;
% perm = reshape(permSmallInit', [1,G.cells.num]);
% rock.perm = perm'*100*milli*darcy;
% rock.perm = perm';

modelHMSmall

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

obj    = objectiveFunction(G, S, W, rock, fluid, simRes);
Jinit  = obj.val;
                         
% adjoint run
adjRes = runAdjoint(simRes, G, S, W, rock, fluid, schedule, controls, ...
                    objectiveFunction, 'VerboseLevel', verboseLevel);
                
grad   = computeGradientDK(simRes, adjRes, G, S, W, rock, fluid, schedule, controls, objectiveFunction);

adjGrad = -cell2mat(grad);

JnormJscale = 2e1 * (adjGrad ./ norm(adjGrad,2));
% JnormJscale = 2 * (adjGrad ./ norm(adjGrad,inf));

% Update parameter
% JnormJscale = abs(JnormJscale);
% rock.perm   = rock.perm + JnormJscale*100*milli*darcy;
% rock.perm   = rock.perm + JnormJscale;

m         = reshape(m, [1,G.cells.num]);
m         = m(:);
m         = m + abs(JnormJscale);
load Kreal;
Kreal     = reshape(Kreal, [G.cells.num,1]);
rock.perm = (Kreal ./ m)*100*milli*darcy;

modelHMSmall

% forward run
simRes = runSchedule(resSolInit, G, S, W, rock, fluid, schedule, ...
                             'VerboseLevel', verboseLevel);

obj    = objectiveFunction(G, S, W, rock, fluid, simRes);
Jpert  = obj.val;
                              
% test to ONE, psi should be ONE if correct !
psi = ( Jpert - Jinit ) ./ (JnormJscale' * adjGrad);  

return;                              
                              
 