% simple2DMSSBSourceExample -- A simple test case comparison between
%                          fine-scale and multiscale Stokes-Brinkman
%                          solvers driven by sources/sinks.  
%
% SYNOPSIS:
%   simple2DMSSBSourceExample
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
physDims  = [0 nx 0 ny];
coarseDim = [5 5];

% generate grid and DOFs
G          = cartGrid2D(cartDims,physDims([2 4]));
G.physDims = physDims;
G          = computeGeometry(G);
Dofs       = findCartDofs(G);

% generate coarse grid
part  = partitionUI(G, coarseDim);
CG    = generateCoarseGrid(G, part);

% generate rock and fluid structs
fluid        = initSingleFluid;
fluid.mu_eff = fluid.mu;
rock.perm    = ones(G.cells.num,1)*darcy()/1000;

% set sources/sinks
csource = G.cells.num/2-ny/2; 
csinks = [1,nx,nx*(ny-1)+1,G.cells.num];

src = addSource([],  csource, ones(size(csource))/numel(csource));
src = addSource(src, csinks, -ones(size(csinks))/numel(csinks));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
Ssb = makeSystemSB(G, Dofs, rock, fluid);

% solve system
[Ssb, FSsb] = solveSystemSB(Ssb, G, Dofs, 'src', src);
FSsb        = nodeToCellData(FSsb, G, Dofs, 'CG', CG, 'part', part);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
CS = generateCoarseSystemSB(G, rock, Dofs, CG, fluid, Ssb.element,...
                             Ssb.basis, 'src', src);

% solve system
[MSsb, CS]    = solveCoarseSystemSB(G, Ssb, CG, CS, Dofs, rock, fluid, part, ...
                                    'src', src);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MSvtot = sqrt(MSsb.v1.^2+MSsb.v2.^2); 
FSvtot = sqrt(FSsb.v1.^2+FSsb.v2.^2); 

integral = [1;4;1;4;16;4;1;4;1];
cvelFS   = sqrt(FSsb.v1(Dofs.Vdofs)*integral/6/6).^2+...
    (FSsb.v2(Dofs.Vdofs)*integral/6/6).^2;
cvelMS   = sqrt(MSsb.v1(Dofs.Vdofs)*integral/6/6).^2+...
    (MSsb.v2(Dofs.Vdofs)*integral/6/6).^2;

fprintf(1,'\nL2 pressure norm: %4.2d', ...  
        norm(abs(FSsb.cellPressure-MSsb.cellPressure))/norm(FSsb.cellPressure));
fprintf(1,'\nL2 velocity norm: %4.2d\n', norm(abs(FSvtot-MSvtot))/norm(FSvtot));

figure('Position',[0 300 800 800]);
subplot(2,2,1); plotCellData(G, cvelFS*day()); cx = caxis;
title('Darcy flux [m/day]'); shading flat; axis equal tight; 
subplot(2,2,2); plotCellData(G, cvelMS*day()); caxis(cx);
title('Stokes-Brinkman flux [m/day]'); shading flat; axis equal tight;
subplot(2,2,3); plotCellData(G, convertTo(FSsb.cellPressure, barsa()));
title('Darcy pressure [bar]'); shading flat; axis equal tight; cx = caxis;
subplot(2,2,4); plotCellData(G, convertTo(MSsb.cellPressure, barsa())); caxis(cx);
title('Stokes-Brinkman pressure [bar]'); shading flat; axis equal tight;
