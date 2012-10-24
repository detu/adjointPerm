%----------------------------------------------------------------------------------
% SYNOPSIS:
%
% DESCRIPTION:
%   Re-Initialize all models after an update
%
% PARAMETERS:
%
% RETURNS:
%   dynamic state and structure arrays for the truth and for all ensemble models
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
       % structural update
       if t == 0 && (iter > 1 || iu > 1)
          if length(nstat) > 4
            load(EnkfFile,'U');
            G = U(sum(nstat(1:4))+1:sum(nstat),:);
          end
       end

       % set up matlab structures for truth (l=1) and ensemble members (l=2)
       for l = 1 : 2

         if l == 1
           ne = 0;
         else
           ne = nmembers;
         end

         for j = min(1,ne) : ne
             
           % construct the model
           if t == 0 || ~isempty(G)
            [grid, well, fluid, rock, p0, s0, phi, K, G] = getModel(workdir,modelname,nmembers,j,G);
           end

           % assign rock properties
           if j > 0
             if iter == 1 && iu == 1
               rock.perm = [K(:,j) K(:,j) K(:,j)/10];
               rock.poro = phi(:,j);
             else
               load(EnkfFile,'U');
               rock.perm = [U(1:nstat(1),j) U(1:nstat(1),j) U(1:nstat(1),j)/10];
               rock.poro = U(nstat(1)+1:sum(nstat(1:2)),j);
               rSol.pressure = U(sum(nstat(1:2))+1:sum(nstat(1:3)),j);
               rSol.s = U(sum(nstat(1:3))+1:sum(nstat(1:4)),j);
               clear U
             end
           end

           if t == 0     
             % initialize solution structures and assemble linear hybrid system
             plss = computeMimeticIP(grid, rock, 'Verbose', true);
             rSol = initResSol (grid, p0, s0);
             rSol.wellSol = initWellSol(well, p0);
             % solve linear system to obtain solution for flow and pressure in reservoir and wells
             rSol = solveIncompFlow(rSol, grid, plss, fluid, 'wells', well);
           end

           if j > 0
             rockE{j} = rock; 
             rSolE{j} = rSol;
             if ~isempty(G)
                 gridE{j} = grid;
             end
             if t == 0, plssE{j} = plss; end
           else
             rockT = rock; 
             rSolT = rSol;
             gridT = grid;
             plssT = plss;
           end

         end

       end

       rock = rockT; 
       rSol = rSolT;
       grid = gridT;
       plss = plssT;
       clear rockT rSolT gridT plssT ne
       