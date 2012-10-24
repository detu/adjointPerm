%----------------------------------------------------------------------------------
% SYNOPSIS:
%
% DESCRIPTION:
%   Simulate synthetic truth and model ensemble to the next update time or to the 
%   specified end time and store production time series.
%
% PARAMETERS:
%
% RETURNS:
%   rSol, rSolE         -   MRST states for truth and ensemble
%   wellData, wellDataE -   production time series
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%----------------------------------------------------------------------------------        

       istep = 1;
       while t < tu   
           
         step = dt; if t + dt > tu, step = tu-t; end
         
         % simulate truth          
         rSol = implicitTransport(rSol, grid, step, rock, fluid, 'wells', well);
         assert(max(rSol.s) < 1+eps && min(rSol.s) > -eps);
         rSol = solveIncompFlow(rSol, grid, plss, fluid, 'wells', well);

         % simulate ensemble
         for j = 1 : nmembers
           if ~isempty(G)
              grd = gridE{j};
           else
              grd = grid;
           end
           rSolE{j} = implicitTransport(rSolE{j}, grd, step, rockE{j}, fluid, 'wells', well);
           assert(max(rSolE{j}.s) < 1+eps && min(rSolE{j}.s) > -eps);
           rSolE{j} = solveIncompFlow(rSolE{j}, grd, plssE{j}, fluid, 'wells', well);
         end
         
         t = t + step;  time = convertTo(t,day);
         
         %---------------------------------------------------------------------------
         % store production data time series (rates in m^3/s, pressures in Pa)
         %---------------------------------------------------------------------------             
         if ~exist('wellData','var')
             for i = 1 : length(rSol.wellSol)
                 wellData(i).loc  = well(i,1).cells;                 
                 wellData(i).bhp  = [];
                 wellData(i).flx  = [];
                 wellData(i).ffl  = [];
                 wellData(i).icv  = [];
                 wellDataE(i).bhp = [];
                 wellDataE(i).flx = [];
                 wellDataE(i).ffl = [];
                 wellDataE(i).icv = [];
             end
             timeSteps = [];
         end       
         for i = 1 : length(rSol.wellSol)
             wellData(i).bhp = [wellData(i).bhp; rSol.wellSol(i).pressure];
             rate = 0;
             for k=1:length(rSol.wellSol(i).flux)
                wellData(i).icv = [wellData(i).icv; rSol.wellSol(i).flux(k)]; 
                rate = rate + rSol.wellSol(i).flux(k);
             end
             wellData(i).flx = [wellData(i).flx; rate];
             
             % compute fractional flow assuming horizontal flow and no 
             % capillary pressure: fw = qw/qt = 1/(1+(muw/krw)*(kro/muo))
             s = rSol.s(well(i).cells);
             kr = fluid.relperm(s); % [krw kro]
             mu = fluid.properties(s); % [muw muo]
             mo = bsxfun(@rdivide, kr, mu); % [krw/muw kro/muo]
             fw = 1/(1 + mo(2)/mo(1));
             if isnan(fw), fw = 0; end
             wellData(i).ffl = [wellData(i).ffl; fw];
         
             if nmembers > 0
               pressure = zeros(1,nmembers);
               fracflow = zeros(1,nmembers);
               flow = zeros(length(rSol.wellSol(i).flux),nmembers);
               for j = 1 : nmembers
                 pressure(1,j) = rSolE{j}.wellSol(i).pressure;
                 s = rSolE{j}.s(well(i).cells);
                 kr = fluid.relperm(s);
                 mu = fluid.properties(s);
                 mo = bsxfun(@rdivide, kr, mu);
                 fw = 1/(1 + mo(2)/mo(1));
                 if isnan(fw), fw = 0; end
                 fracflow(1,j) = fw;
               end
               wellDataE(i).bhp = [wellDataE(i).bhp; pressure];
               for k = 1 : length(rSol.wellSol(i).flux)
                 for j = 1 : nmembers
                    flow(k,j) = rSolE{j}.wellSol(i).flux(k);
                 end
                 wellDataE(i).icv = [wellDataE(i).icv; flow(k,:)];
               end
               wellDataE(i).flx = [wellDataE(i).flx; sum(flow,1)];
               wellDataE(i).ffl = [wellDataE(i).ffl; fracflow];
             end
         end
         
         if t > 0
            rSol1 = rSol;
         end         
         
         timeSteps = [timeSteps; time]; istep = istep + 1;
         clear rate flow pressure fracflow
            
       end
       
       disp(['finished all simulations. time = ' num2str(time)])      
       