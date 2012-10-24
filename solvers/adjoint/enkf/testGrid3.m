%% Real field model example
%
% This example uses the Eclipse grid for the Norne field made available by
% Statoil. Random permeability realizations are generated.
%

%linsolve = @(S,b) agmg(S, b, [], 1e-10, [], 1);
linsolve = @mldivide;


%% ONLY FOR TESTING
% ne = 1; j = ne;

%% Check for existence of input model data
grdecl = fullfile(ROOTDIR, 'examples\data\','NORNE','NORNE.GRDECL');

if ~exist(grdecl, 'file'),
   error('Model data is not available.')
end


%% Read and process the model
grdecl = readGRDECL(grdecl)

% select first layer only
actnum = reshape(grdecl.ACTNUM,grdecl.cartDims); actnum(:,:,2:end) = 0; %actnum(:,:,11:end) = 0;
grdecl.ACTNUM = reshape(actnum,prod(grdecl.cartDims),1); clear actnum

G = processGRDECL(grdecl); 
G = computeGeometry(G(1));

faces = find(G.faces.tag>0);
actnum = grdecl.ACTNUM; clear grdecl;

%% Set rock and fluid data

gravity off

for jj = 1 : ne + 1
    kx = logNormLayers(G.cartDims, 100*[1:G.cartDims(3)], 'sz', [11, 11, 1]);

    kx = kx(G.cells.indexMap);
    kx = 2 + (kx-min(kx))/(max(kx)-min(kx))*4000;
   
    K(:,jj) = kx;
    phi(:,jj) = 0.25*ones(size(kx));
    
end

K = convertFrom(K, milli*darcy);

rock.perm = [K(:,end) K(:,end) K(:,end)/10];
rock.poro = phi(:,end);
                         
fluid = initCoreyFluid('mu' , [   0.4,  0.9]*centi*poise, ...
                         'rho', [1014, 859]*kilogram/meter^3, ...
                         'n'  , [   4,   4]                 , ...
                         'sr' , [ 0.2, 0.2]                 , ...
                         'kwm', [   1,   1]);                         

   
%% Initial state
p0 = 350*barsa(); %400*barsa();
s0 = 0.2;


%% Introduce wells

% Set vertical injectors
nz = G.cartDims(3);
I = [ 9, 26,  8, 25, 35, 10];
J = [14, 14, 35, 35, 68, 75];
% R = 0.25 * [ 4,  4,  4,  4,  4,  4]*1000*meter^3/day;
bhp = 350*barsa();
nIW = 1:numel(I); W = [];
for i = 1 : numel(I),
%    W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
%                     'Val', R(i), 'Radius', 0.1, 'Comp_i', [1,0,0], ...
%                     'name', ['I$_{', int2str(i), '}$']);
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
                    'Val', p0, 'Radius', 0.1, 'Comp_i', [1,0,0], ...
                    'name', ['I$_{', int2str(i), '}$']);                
end

% Set vertical producers
I = [17, 12, 25, 35, 15];
J = [23, 51, 51, 95, 94];
bhp = 335*barsa();
nPW = (1:numel(I))+max(nIW);
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
                    'Val', bhp, 'Radius', 0.1, ...
                    'name', ['P$_{', int2str(i), '}$']);
end

%% Initialize solution structures and assemble linear hybrid system
% plss = computeMimeticIP(G, rock, 'Verbose', true); % check if this is different for modified structures, otherwise use one only
% rSol = initResSol (G, p0, s0);
% rSol.wellSol = initWellSol(W, bhp);

%% Solve linear system to obtain solution for flow and pressure in reservoir and wells
% rSol = solveIncompFlow(rSol, G, plss, fluid, 'wells', W); %,'LinSolve',linsolve);

%% Simulate
% Tend           = 15.00*year();      % end time
% dT             = 0.50*year();       % simulation timestep size 0.01*year() (30*day())
%  t = 0; time = convertTo(t,day);
% while t < Tend   
% 
%     step = dT; if t + dT > Tend, step = Tend-t; end
% 
%     % simulate truth          
%     rSol = implicitTransport(rSol, G, step, rock, fluid, 'wells', W);
%     assert(max(rSol.s) < 1+eps && min(rSol.s) > -eps);
%     rSol = solveIncompFlow(rSol, G, plss, fluid, 'wells', W);
%     
%     t = t + step; time = convertTo(t,day);
%     
% end


%% ONLY FOR TESTING
% figure(j)
% axis equal
% % plotCellData(G,convertTo(rock.perm(:,1),milli*darcy),'EdgeColor','k','FaceAlpha', 1.0);
% % plotCellData(G,convertTo(K(:,j),milli*darcy),'EdgeColor','k','FaceAlpha', 1.0);
% plotCellData(G,rSol.s,'EdgeColor','k','FaceAlpha', 1.0);
% plotWell(G, W, 'height', 75, 'color', 'k');
% axis off, view(15,60), h=colorbar('horiz');
% % caxis([0 2000]);
% zoom(2.5), title('permeability [mD]');
% axis tight off, view([-40, 30]);


