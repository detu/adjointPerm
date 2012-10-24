%----------------------------------------------------------------------------------
% SYNOPSIS:
%   mrstEnKF
%
% DESCRIPTION:
%   Basic Ensemble Kalman Filter framework for the SINTEF reservoir simulator MRST.
%   The code has been tested with Matlab version R2011a.
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
% This version: 2012a.
%---------------------------------------------------------------------------------- 
 cdd = pwd;
 addpath(genpath(cdd));

%---------------------------------------------------------------------------------- 
% input section
%----------------------------------------------------------------------------------
%  inputSettings;
inputSettings1;

%---------------------------------------------------------------------------------- 
% loop over experiments
%----------------------------------------------------------------------------------           
 for ir = 1 : nrepeats + 1

 rng(ir,'v5normal'); rng(ir,'v5uniform');     
 
 %--------------------------------------------------------------------------------- 
 % preparation
 %---------------------------------------------------------------------------------
%  if exist(workdir,'dir')
%    [~, ~, ~] = rmdir(workdir,'s');
%  end
%  [status, message, messageid] = mkdir(workdir);
%  cd(workdir);

 %---------------------------------------------------------------------------------
 % loop over update times
 %---------------------------------------------------------------------------------
 t = 0; 
 time = convertTo(t,day);
 tda = union(dtu:dtu:tend,tend);
 tdo = union(dto:dto:tend,tend);
 tda = union(tda,tsim);
 for iu = 1 : length(tda)
   tu = tda(iu);
   
   % select observation times for the current update step
   if iu > 1,
      obsTimes = intersect(tdo(tdo > tda(iu-1)),tdo(tdo <= tda(iu)));
   else
      obsTimes = tdo(tdo <= tda(iu));
   end
   obsTimes = convertTo(obsTimes,day);    
   
   %-------------------------------------------------------------------------------
   % loop over iterations  
   %-------------------------------------------------------------------------------
   iter = 1; niter = niterations + 1;
   if (restart == 0 && niterations == 1)
       niter = niterations;
   end
   while iter <= niter

     % set the start time for simulation
     if iter > 1
         t = tda(max(1,iu-1));
         clear timeSteps wellData*
     end
     if restart == 1
         t = 0;
         clear timeSteps wellData*
     end
     time = convertTo(t,day);

     %-----------------------------------------------------------------------------
     % initialize the model (if t = 0 -> grid, well, fluid, rock, structure
     %------------------------------------------------------------------------------
     initializeModel;

     %------------------------------------------------------------------------------
     % simulate to next update time
     %------------------------------------------------------------------------------
     simulateModel;
     
     %------------------------------------------------------------------------------
     % construct ensemble state matrix A
     %------------------------------------------------------------------------------
     prm = rock.perm(:,1);
     por = rock.poro;
     sat = rSol.s;
     prf = rSol.pressure;
     xt = [prm; por; prf; sat];
     if nmembers > 0
         for j=1:nmembers
             prme(:,j) = rockE{j}.perm(:,1);
             pore(:,j) = rockE{j}.poro;
             prfe(:,j) = rSolE{j}.pressure;
             sate(:,j) = rSolE{j}.s;
         end
         A = [prme; pore; prfe; sate];
     end 
     nstat = [length(prm),length(por),length(prf),length(sat)];
     clear prm* por* prf*
     
     if ~isempty(G)
         A = [A; G];
         nstat = [nstat, size(G,1)];
     end
     
     %------------------------------------------------------------------------------
     % observations
     %------------------------------------------------------------------------------
     if isempty(find(strncmp('sat',obstype,3),1))
        [Y,R,D,loc] = makeObservations(wellData,wellDataE,timeSteps,obsTimes,...
                      obstype,transformobs,sigma);
     else
        [Y,R,D,loc] = makeObservations(wellData,wellDataE,timeSteps,obsTimes,...
                      obstype,transformobs,sigma,sat,sate);  
     end
     clear sat*
     
     nobs = size(Y,1);
     if nobs == 0, iter = niter; end
     
     %------------------------------------------------------------------------------
     % model-data misfit
     %------------------------------------------------------------------------------
     if nobs > 0 && nmembers > 0
       E = D - Y; 
       mse = 0;
       for i = 1 : nmembers
         for j = 1 : size(E,1)  
           mse = mse + E(j,i) * E(j,i) / R(j);
         end
       end
       mse = mse / nmembers;
     end
     clear E
     
     %------------------------------------------------------------------------------
     % update the models if needed
     %------------------------------------------------------------------------------
     if nmembers > 0 && iter <= niterations
         
       %----------------------------------------------------------------------------
       % transformation of states
       %----------------------------------------------------------------------------  
       A = applyTransform(A,nstat,transform);
       
       %----------------------------------------------------------------------------
       % initialize iterative updating schemes
       %----------------------------------------------------------------------------
       if iter == 1
         A0 = A; J = []; iterations = []; beta = 1.0; Jmin = 1.e+32;
         EnkfFile = ['states_' num2str(time) '_1.mat'];        
         if exist(EnkfFile,'file'), delete(EnkfFile); end       
       else
         EnkfFile = ['states_' num2str(time) '_' num2str(iter-1) '.mat']; 
         load(EnkfFile,'A0','A1','J','beta'); % 'Y',
         Jmin=min(J);
         EnkfFile = ['states_' num2str(time) '_' num2str(iter) '.mat'];
       end
                            
       epsilon = [1.e-6 1.e-6];                        
       if (mse <= size(Y,1)) || (iter > 1 && ( mse < Jmin && Jmin-mse < epsilon(2)*Jmin))
         iter = niter;
       end
       if mse > Jmin
         U = A1;
         beta = 0.5*beta;
         if beta < 0.05
           iter = niter;
         end
       else
         J = [J mse];
         A1 = A; 
         beta = min(1.25*beta, 1);
                
         if strcmp(method,'EnRML') 
                
           U = EnRML(D, Y, R, A, A0, beta);  
                      
         elseif strcmp(method,'EnKF')
             
           U = EnKF(D, Y, R, A, beta);  
         
         elseif strcmp(method,'LEnKF')
             
           LEnKF;
           
         else  
                   
           error(' *** method not yet implemented ***');
           
         end    

       end
       
       %----------------------------------------------------------------------------
       % inflation
       %----------------------------------------------------------------------------
       if inflation > 0
         inflateEnsemble;
       end
       
       %----------------------------------------------------------------------------
       % back transformation of states
       %----------------------------------------------------------------------------
       U = applyTransform(U, nstat, itransform);
           
       %----------------------------------------------------------------------------
       % fix updated states
       %----------------------------------------------------------------------------
       U = fixStates(U,nstat);

       %----------------------------------------------------------------------------
       % store iteration results
       %----------------------------------------------------------------------------
       save(EnkfFile,'xt','U','A0','Y','A1','J','beta','grid');
       
     else
         
       if nmembers > 0
         EnkfFile = ['states_' num2str(time) '_' num2str(iter) '.mat']; 
         U = A;
         if isempty(G)
            save(EnkfFile,'xt','U','grid');
         else
            save(EnkfFile,'xt','U','grid','gridE');
         end
         clear U
       end
        
     end
        
     %------------------------------------------------------------------------------
     % save well data
     %------------------------------------------------------------------------------
     summary = ['summary_' num2str(time) '_iter' num2str(iter) '.mat'];
     if nmembers > 0
       save(summary,'wellData','wellDataE','timeSteps','settings');
     else
       save(summary,'wellData','timeSteps','settings');
     end
     
     iter = iter+1;
     
   end
 
 end
 

 cd(cdd);
 
%----------------------------------------------------------------------------------
 end % repeated experiments
%----------------------------------------------------------------------------------