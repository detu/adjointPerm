% simple3DSBSourceExample -- A simple test case of Stokes-Brinkman solvers
%                            driven by sinks/sources. 
%
% SYNOPSIS:
%   simple3DSBSourceExample
%
% PARAMETERS:
%   None
%
% RETURNS:
%   Nothing, though being implemented as a script, the variables remain
%   available in the base workspace upon completion.
 
% set grid parameters
nx        = 9;  
ny        = 9;
nz        = 4;
cartDims  = [nx ny nz];
physDims  = [0 nx 0 ny 0 nz];

% generate grid and DOFs
G          = cartGrid(cartDims,physDims([2 4 6]));
G.physDims = physDims;
G          = computeGeometry(G);
Dofs       = findCartDofs(G);

% generate rock and fluid structs
fluid        = initSingleFluid;
fluid.mu_eff = fluid.mu;;
rock.perm    = ones(G.cells.num,1)*darcy()/1000;

% set sources/sinks
csource = [1, nx*(ny-1)+1, nx*ny*nz-nx+1, nx*ny*(nz-1)+1]; 
csinks = [nx, nx*ny, nx*ny*(nz-1)+nx, G.cells.num];

src = addSource([],  csource, ones(size(csource))/numel(csource));
src = addSource(src, csinks, -ones(size(csinks))/numel(csinks));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
Ssb = makeSystemSB(G, Dofs, rock, fluid);

% solve system
[Ssb, FSsb] = solveSystemSB(Ssb, G, Dofs, 'src', src);

p = mean(FSsb.p(Dofs.Pdofs),2);

integral=[1;4;1;4;16;4;1;4;1;4;16;4;16;64;16;4;16;4;1;4;1;4;16;4;1;4;1];
vxc      = FSsb.v1(Dofs.Vdofs)*integral/6/6/6; 
vyc      = FSsb.v2(Dofs.Vdofs)*integral/6/6/6; 
vzc      = FSsb.v3(Dofs.Vdofs)*integral/6/6/6; 
vtot     = sqrt(vxc.^2+vyc.^2+vzc.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[0 300 800 400]);
subplot(1,2,1); plotCellData(G, log10(vtot*day())); view(3)
title('Log10 Velocity [m/day]'); shading flat; axis equal tight;
subplot(1,2,2); plotCellData(G, log10(convertTo(p, barsa()))); view(3)
title('Log10 Pressure [bar]'); shading flat; axis equal tight;

