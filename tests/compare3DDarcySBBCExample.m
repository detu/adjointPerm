% compare2DDarcySBBCExample -- A simple test case comparison between
%                              Darcy and Stokes-Brinkman solvers driven
%                              by pressure BCs.
%
% SYNOPSIS:
%   compare2DDarcySBBCExample
%
% PARAMETERS:
%   None
%
% RETURNS:
%   Nothing, though being implemented as a script, the variables remain
%   available in the base workspace upon completion.

% set grid parameters
nx        = 5;  
ny        = 5;
nz        = 20;
cartDims  = [nx ny nz];
physDims  = [0 nx 0 ny 0 nz];

% generate grid and DOFs
G          = cartGrid(cartDims,physDims([2 4 6]));
G.physDims = physDims;
G          = computeGeometry(G);
Dofs       = findCartDofs(G);

% generate rock and fluid structs
fluid        = initSingleFluid;
fluid.mu_eff = 0;
rock.perm    = ones(G.cells.num,1)*darcy()/1000;
rock.poros   = repmat(0.3, [G.cells.num, 1]);
gravity off;

% set BCs
ind         = any(G.faces.neighbors==0,2);
bcfaces     = find(ind);
dum         = false(G.faces.num,1); 
dum(bcfaces)= true;
tags        = G.cellFaces(dum(G.cellFaces(:,1)),2);
bcfaces     = G.cellFaces(dum(G.cellFaces(:,1)),1);

face_le = bcfaces(tags==1); face_r = bcfaces(tags==2);
face_l  = bcfaces(tags==3); face_u = bcfaces(tags==4);
face_b  = bcfaces(tags==5); face_t = bcfaces(tags==6);
bcfaces = [face_l; face_u; face_le; face_r];

p_t = 2*barsa(); p_b = 1*barsa();

BCsb = addBCSB([],   face_t, 'pressure',   repmat(p_t, numel(face_t), 1), ...
               G, Dofs);
BCsb = addBCSB(BCsb, face_b , 'pressure',   repmat(p_b,  numel(face_b),  1), ...
               G, Dofs);
BCsb = addBCSB(BCsb, bcfaces, 'velocity_n', repmat(0,    numel(bcfaces), 1), ...
               G, Dofs); 
BCd  = addBC(  [],   face_t, 'pressure',   repmat(p_t, numel(face_t), 1));
BCd  = addBC(  BCd,  face_b , 'pressure',   repmat(p_b, numel(face_b),  1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Stokes-Brinkman (Taylor-Hood)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
Ssb = makeSystemSB(G, Dofs, rock, fluid);

% solve system
[Ssb, FSsb] = solveSystemSB(Ssb, G, Dofs, 'bc', BCsb);
FSsb        = nodeToCellData(FSsb, G, Dofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fine-scale Darcy (Raviart-Thomas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble system
S = computeMimeticIP(G, rock);

% solve system
[FSd, xwRef] = solveIncompFlow(initResSol(G, 0.0), initWellSol([], 0.0), ...
                                 G, S, fluid, 'bc', BCd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
fluxFSd  = accumarray(cellNo, abs(FSd.cellFlux));
fluxFSsb = accumarray(cellNo, abs(FSsb.cellFlux));

fprintf(1,'\nL2 pressure norm: \t%4.2d', ...
        norm(abs(FSsb.cellPressure-FSd.cellPressure))/norm(FSd.cellPressure));

fprintf(1,'\nL2 flux norm:\t\t%4.2d\n',...
        norm(abs(fluxFSsb-fluxFSd))/norm(fluxFSd));

figure('Position',[0 300 800 800]);
subplot(2,2,1); plotCellData(G, log10(fluxFSd*day()));cx = caxis; view(3);
title('Log10 Darcy flux [m/day]'); shading flat; axis equal tight; 
subplot(2,2,2); plotCellData(G, log10(fluxFSsb*day())); caxis(cx); view(3);
title('Log10 Stokes-Brinkman flux [m/day]'); shading flat; axis equal tight;
subplot(2,2,3); plotCellData(G, convertTo(FSd.cellPressure, barsa()));
view(3);
title('Darcy pressure [bar]'); shading flat; axis equal tight; cx = caxis;
subplot(2,2,4); plotCellData(G, convertTo(FSsb.cellPressure, barsa()));
caxis(cx); view(3);
title('Stokes-Brinkman pressure [bar]'); shading flat; axis equal tight; 