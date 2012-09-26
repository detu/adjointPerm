function plot_num = findPlotNum(Nx, Ny, Vdofs)
% findPlotNum - Finds the correct order of elements for plotting the
%               velocity defined by the Taylor-Hood elements
%
% SYNOPSIS:
%   plot_num = findPlotNum(Nx, Ny, Vdofs)
%
% PARAMETERS:
%   Nx    - The number of cells in x-direction
%
%   Ny    - The number of cells in y-direction
%
%   Vdofs - A matrix with the DOFs for the velocities listed per cell (Dofs.Vdofs) 
%
% RETURNS:
%   plot_num - The correct order of elements for plotting the SB velocity

cell_num=Nx*Ny;

plot_num=zeros(2*Ny+1,2*Nx+1);
i=1:2:2*Ny;
j=1:2:2*Nx;

plot_num(i,j)=reshape(Vdofs(:,1),length(j),length(i))';
plot_num(i,Nx*2+1)=Vdofs(Nx:Nx:cell_num,3);
plot_num(Ny*2+1,j)=Vdofs(Ny:Ny:cell_num,3);
plot_num(Ny*2+1,j)=Vdofs(cell_num-Nx+1:1: ...
                                         cell_num,7);
plot_num(Ny*2+1,2*Nx+1)=Vdofs(cell_num,9);

plot_num(i+1,j)=reshape(Vdofs(:,4),length(j),length(i+1))';
plot_num(i+1,end)=Vdofs(Nx:Nx:cell_num,6);

plot_num(i,j+1)=reshape(Vdofs(:,2),length(j),length(i+1))';
plot_num(end,j+1)=Vdofs(cell_num-Nx+1:cell_num,8);
plot_num(i+1,j+1)=reshape(Vdofs(:,5),length(j+1),length(i+1))';