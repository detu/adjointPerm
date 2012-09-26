% simple2DSBBCExample -- A simple test case of Stokes-Brinkman solvers
%                        driven by pressure BCs. 
%
% SYNOPSIS:
%   simple2DSBBCExample
%
% PARAMETERS:
%   None
%
% RETURNS:
%   Nothing, though being implemented as a script, the variables remain
%   available in the base workspace upon completion.

% set grid parameters
nx        = 30;    
ny        = 30;  
cartDims  = [nx ny];
physDims  = [0 1 0 1];

% generate grid and DOFs
G          = cartGrid2D(cartDims,physDims([2 4]));
G.physDims = physDims;
G          = computeGeometry(G);
Dofs       = findCartDofs(G);

% generate rock and fluid structs
fluid        = initSingleFluid;
fluid.mu_eff = fluid.mu;;
rock.perm    = ones(G.cells.num,1)*darcy()/1000;

% set BCs
ind         = any(G.faces.neighbors==0,2);
bcfaces     = find(ind);
dum         = false(G.faces.num,1); 
dum(bcfaces)= true;
tags        = G.cellFaces(dum(G.cellFaces(:,1)),2);
bcfaces     = G.cellFaces(dum(G.cellFaces(:,1)),1);

face_le = bcfaces(tags==1); face_r = bcfaces(tags==2);
face_l  = bcfaces(tags==3); face_u = bcfaces(tags==4);
bcfaces = [face_l;face_u];

p_le = 2*barsa(); p_r = 1*barsa();

BCsb = addBCSB([],   face_le, 'pressure',   repmat(p_le, numel(face_le), 1), G, Dofs);
BCsb = addBCSB(BCsb, face_r , 'pressure',   repmat(p_r,  numel(face_r),  1), G, Dofs);
BCsb = addBCSB(BCsb, bcfaces, 'velocity_n', repmat(0,    numel(bcfaces), 1), G, Dofs); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
Ssb = makeSystemSB(G, Dofs, rock, fluid);

% solve system
[Ssb, FSsb] = solveSystemSB(Ssb, G, Dofs, 'bc', BCsb);

p = mean(FSsb.p(Dofs.Pdofs),2);

integral = [1;4;1;4;16;4;1;4;1];
vxc      = FSsb.v1(Dofs.Vdofs)*integral/6/6; 
vyc      = FSsb.v2(Dofs.Vdofs)*integral/6/6; 
vtot     = sqrt(vxc.^2+vyc.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[0 300 800 400]);
subplot(1,2,1); plotCellData(G, vtot*day());
title('Velocity [m/day]'); shading flat; axis equal tight;
subplot(1,2,2); plotCellData(G, convertTo(p, barsa()));
title('Pressure [bar]'); shading flat; axis equal tight;

%figure('Position',[0 300 400 400]);
%plot_streamlines(G, Dofs, FSsb.v1, FSsb.v2, 30, 'direction', 'all'); 
%title('Log10 (total vel per cell) + streamlines');
