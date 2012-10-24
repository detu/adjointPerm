
%% input
workDir = [ROOTDIR 'TEST3_HM\'];
model   = 'testGrid3';
time    = 16425;
iter    = 2;

%% load grid
load([workDir model '.mat']);
ng   = size(K,1);
ns   = numel(G);
nstat = [ng ng ng ng]; if ns > 0, nstat = [nstat ns]; end

% axis equal
% plotCellData(grid, convertTo(rock.perm(:,1),milli*darcy), 'EdgeColor', 'k', 'FaceAlpha', 1.0);
% plotWell(grid, well, 'height', 75, 'color', 'k');
% axis tight off, view(15,60), h=colorbar('horiz');
% zoom(2.5), title('permeability [mD]');

%% load experiment settings
load([workDir 'summary_' num2str(time) '_iter' num2str(iter) '.mat']);

for i = 1: length(settings)
    itransform{i} = settings(i).transform;
end
for i = 1 : length(itransform)
     itransform{i}(1) = -1 * itransform{i}(1);
end
ne   = settings.nmembers;

%% construct states
load([workDir 'states_' num2str(time) '_' num2str(iter) '.mat']);

% update states
grd = U;

% initial states (may still need to be transformed back)
% grd = transformStates(A0,nstat,itransform);

prm = xt(1:ng); por = xt(ng+1:2*ng); prf = xt(2*ng+1:3*ng); sat = xt(3*ng+1:4*ng);
for j = 1 : ne
    prme(:,j) = grd(1:ng,j); 
    pore(:,j) = grd(ng+1:2*ng,j); 
    prfe(:,j) = grd(2*ng+1:3*ng,j); 
    sate(:,j) = grd(3*ng+1:4*ng,j);
end

%% plotting
nr = 2; nc = 2; em = [1 2 3];
subplot(nr,nc,1)
axis equal
plotCellData(grid, sat, 'EdgeColor', 'k', 'FaceAlpha', 1.0);
plotWell(grid, well, 'height', 75, 'color', 'k');
axis tight off, view(15,60), h=colorbar('horiz'); zoom(2)
for j = 1:length(em)
    subplot(nr,nc,j+1)
    axis equal
    plotCellData(grid, sate(:,em(j)), 'EdgeColor', 'k', 'FaceAlpha', 1.0);
    plotWell(grid, well, 'height', 75, 'color', 'k');
    axis tight off, view(15,60), h=colorbar('horiz'); zoom(2)
end


%% printing

% set(gcf,'PaperUnits','normalized')
% set(gcf,'PaperType','a4letter')
% set(gcf,'PaperPosition',[.01 .01 .99 .99])
% set(gcf,'PaperOrientation','Landscape')
% figname=fullfile(pwd,['plotGrid.pdf']);
% print(gcf,'-dpdf',figname);


