% simple2DMSSBBCExample -- A simple test case comparison between
%                          fine-scale and multiscale Stokes-Brinkman
%                          solvers driven by pressure BCs.  
%
% SYNOPSIS:
%   simple2DMSSBBCExample
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
FSsb        = nodeToCellData(FSsb, G, Dofs, 'CG', CG, 'part', part);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
CS = generateCoarseSystemSB(G, rock, Dofs, CG, fluid, Ssb.element,...
                             Ssb.basis, 'bc', BCsb);

% solve system
[MSsb, CS]    = solveCoarseSystemSB(G, Ssb, CG, CS, Dofs, rock, fluid, part, ...
                                    'bc', BCsb);
 
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
