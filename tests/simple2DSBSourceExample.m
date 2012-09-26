% simple2DSBSourceExample -- A simple test case of Stokes-Brinkman solvers
%                            driven by sinks/sources. 
%
% SYNOPSIS:
%   simple2DSBSourceExample
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

% set sources/sinks
csource = G.cells.num/2-ny/2; 
csinks = [1,nx,nx*(ny-1)+1,G.cells.num];

src = addSource([],  csource, ones(size(csource))/numel(csource)/day());
src = addSource(src, csinks, -ones(size(csinks))/numel(csinks)/day());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
Ssb = makeSystemSB(G, Dofs, rock, fluid);

% solve system
[Ssb, FSsb] = solveSystemSB(Ssb, G, Dofs, 'src', src);

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
