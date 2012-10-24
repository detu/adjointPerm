%%-------------------------------------------------------------------------
% EnKF implementation using MRST
% This has been done by Andrea - DTU on his visit at NTNU January 2012
%--------------------------------------------------------------------------

clear all
close all
clc

Ne = 6; % number of ensemble realizations
trueRealization = 1; % true realization in the ensemble

% initial setup, inital true state xt, initial ensemble items Xap, initial input U,

%we consider equals pressure and saturations for all the ensemble items
[xt, Xap, U, auxdata, R, H] = setupSimulation(Ne,trueRealization);
xam = mean(Xap,2); %mean of the initial ensembles
nassim = 1;   % number of assimilation cycles

nX = size(xt,1);
nY = size(H,1);

% background ensemble states
Xfp = zeros(nX,Ne);

% Root mean square of the full state vector for each ensemble and mean
rms_Xap = zeros(Ne,nassim);
rms_xam = zeros(1,nassim);

%loading parameters
%[G, rock, fluid,wellPar,simParam,modelPar] = deal(auxdata{:});
G = deal(auxdata{1});
numCells = G.cells.num;


%% Initial Statistics
% save mean permeability field 
permMean = zeros(numCells,nassim+1);

% Root mean square of the permeability field for each ensemble and mean
rms_permEn = zeros(Ne,nassim+1);
rms_permMean = zeros(1,nassim+1);

% Root mean square of the full state vector for each ensemble items
for k2 = 1:Ne
    rms_Xap(k2,1) = norm((Xap(:,k2)-xt),2)/sqrt(nX);
end
% Root mean square of the full state vector for the mean estimate
rms_xam(1,1) = norm((xam-xt),2)/sqrt(nX);

% Root mean square of the permeability field for each memeber
for k2 = 1:Ne
    rms_permEn(k2,1) = norm((Xap((nY+1):( nY +numCells),k2)-xt((nY+1):(nY +numCells),1)),2)/sqrt(numCells);
end
% Root mean square of the permeability field
rms_permMean(1,1) = norm((xam((nY+1):(nY +numCells),1)-xt((nY+1):(nY +numCells),1)),2)/sqrt(numCells);

% save mean permeability field
permMean(:,1) = xam((nY+1):(nY +numCells),1);



%% Assimilation
disp(['Cycling ON the reservoir...']);
for k1=1:nassim

  disp(['assimilation step ' int2str(k1)])

  % step to the next assimilation time
  Xb = Xap; 
  
  % advance truth with the full nonlinear model
  [xt outputData] = forwardSimulation(xt, U, auxdata);
  
%   [tims states] = ode45(@derivsL63,[0 ntimes],xt,[],par); 
%   xt = states(end,:)';  % next true state

  
  % Make observation on the truth (without added noise)
  Yt = H * xt; %Yt = my_makeObservation(xt,...)
  
  
  % advance background ensemble with the full nonlinear model
  for k2 = 1:Ne
	  [Xfp(:,k2)] = forwardSimulation(Xb(:,k2), U, auxdata);
  end
  xfm = mean(Xfp,2); %mean of the forecasted states
  
  

  % new background error covariance from ensemble estimate
  % we don't consider errors in the model
  
  
  
  % update step for ensemble mean (xam), assimilated states (Xap), and
  % assimilated covariance error Pa
  [xam, Xap, Pa] = enkfUpdate(xfm, Xfp, Yt, H, R);


  

  
  %% Save statistics 
  
  % Root mean square of the full state vector for each ensemble items
  for k2 = 1:Ne
    rms_Xap(k2,k1+1) = norm((Xap(:,k2)-xt),2)/sqrt(nX);
  end  
  % Root mean square of the full state vector for the mean estimate
  rms_xam(1,k1+1) = norm((xam-xt),2)/sqrt(nX);
 
  % Root mean square of the permeability field for each memeber
  for k2 = 1:Ne
    rms_permEn(k2,k1+1) = norm((Xap((nY+1):( nY +numCells),k2)-xt((nY+1):(nY +numCells),1)),2)/sqrt(numCells);
  end  
  % Root mean square of the permeability field
  rms_permMean(1,k1+1) = norm((xam((nY+1):(nY +numCells),1)-xt((nY+1):(nY +numCells),1)),2)/sqrt(numCells);
  
  % save mean permeability field 
  permMean(:,k1+1) = xam((nY+1):(nY +numCells),1);

      
  
end % end of assimilation loop

[simRes,G] = deal(outputData{:});
figure(91)
clf
subplot(1,3,1)
plotCellData(G,log10(xt((nY+1):(nY +numCells),1)));axis tight; axis square;
title('True')

subplot(1,3,2)
plotCellData(G,log10(permMean(:,1)));axis tight; axis square;
title('Initial Estimated mean')

subplot(1,3,3)
plotCellData(G,log10(permMean(:,end)));axis tight; axis square;
title('Final Estimated mean')

%plotWell(G,W); axis tight; axis square;



a = 1
