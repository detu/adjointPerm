function normS = testThetaScheme(theta)
% twophaseTransport -- Examples of two-phase transport on small Cartesian
%                      grid with homogeneous, isotropic permeability.
%                      Compares transport solvers using explicit, implicit
%                      and theta scheme. 
% SYNOPSIS:
%           testThetaScheme()
%           testThetaScheme(theta)
%   normS = testtehtaScheme(...)
%
% DESCRIPTION:
%   Function testThetaScheme runs the following example:
%     1) Water flooding in a quarter five-spot well configuration on a
%        10-by-10 Cartesian grid with homogeneous, isotropic permeability
%        and porosity.  The reservoir is initially saturated with oil.
%
%     
% PARAMETERS:
%   theta - Paramenter in theta-scheme.   
%
% RETURNS:
%   normS - Euclidian norm of saturation differences, NORM(s_e - s_i), with
%           's_e' being the saturations computed with the explicit
%           transport solver and 's_i' being the saturations computed with
%           the implicit transport solver.  One scalar value for each test.

% $Id: testThetaScheme.m $

if nargin == 0,
   theta = 0.5;    
end

doPlot  = true;
verbose = true;
verbose_solvers = true; % Display output from transport solvers.

if verbose  
   disp('** Running test of two-phase transport **');
   tic
end

%--------------------------------------------------------------------------
%- Initialize system ------------------------------------------------------
%
dims       = [10, 1, 10];
G          = cartGrid(dims);
G          = computeGeometry(G, 'Verbose', verbose);
rock.poro  = repmat(0.3, [G.cells.num, 1]);
rock.perm  = repmat(100, [G.cells.num, 1]) .* darcy() / 1000;
fluid      = initSimpleFluid();

gravity('off');

S = computeMimeticIP(G, rock, 'Verbose', verbose);

normS = zeros([2, 1]);


%- Quarter 5-spot --------------------------------------------
% Test with initial water saturation = 0 and injection of water.
%
if verbose,
   fprintf([ '\n\n****************************\n', ...
                 '** Test 1: Quarter 5-spot **\n', ...
                 '****************************\n']);
   tic
end
sourceCells = [ 1 ; G.cells.num ];
sources     = [ 1 ;    -1       ] ./ 86400;  % m^3/s
src         = addSource([], sourceCells, sources, 'sat', [1,0,0; 1,0,0]);
pv          = poreVolume(G,rock);
tf          = sum(pv) / 2;
timestep    = tf / 5;

t = 0;

tf       = tf       * day(); % Transport solvers expect time in seconds
timestep = timestep * day();

%-----------------------------------------------------------------------
%- Compute initial solution of pressure equation -----------------------
%
resSol = initResSol(G, 0.0, [0,0]);

resSol      = solveIncompFlow(resSol, [], G, S, fluid, 'src', src);
resSol_impl = resSol;
resSol_t    = resSol;

if doPlot,
   % Plot inital saturation
   figure;

   subplot(3,1,1),
   h = plotCellData(G, resSol.s(:,1));
   view([0, 0]);
   axis tight; axis equal;

   subplot(3,1,2),
   h2 = plotCellData(G, resSol.s(:,1));
   view([0, 0]);
   axis tight; axis equal;
   
   
   subplot(3,1,3),
   h3 = plotCellData(G, resSol.s(:,1));
   view([0, 0]);
   axis tight; axis equal;
   
   
end

%-----------------------------------------------------------------------
%- Transport loop - solve saturation equation and update pressures -----
%
while t < tf,
   
   resSol      = explicitTransport(resSol, [], G, timestep, ...
                                   rock, fluid, 'src', src, ...
                                   'verbose', verbose_solvers);

   resSol_impl = implicitTransport(resSol_impl, [], G, timestep, ...
                                   rock, fluid, 'src', src,      ...
                                   'verbose', verbose_solvers);


   resSol_t = thetaTransport(resSol_t, [], G, timestep, ...
                             rock, fluid, 'theta', theta, 'src', src, ...
                             'verbose', verbose_solvers);


   if max(max(abs([resSol.s(:,1), resSol_impl.s(:,1), ...
                  resSol_t.s(:,1)]))) > 1 + 1.0e-5,
      disp('ERROR: ******* Saturation exceeds 1 **********')
      break
   end

   if doPlot,
      % Plot saturation
      delete(h);
      subplot(3,1,1)
      h = plotCellData(G, resSol.s(:,1));
      title(['Test ', int2str(i), ...
         '. Water saturation - explicit solver'])
      view([0, 0]), axis tight equal
      caxis([0, 1]), colorbar

      delete(h2);
      subplot(3,1,2),
      h2 = plotCellData(G, resSol_impl.s(:,1));
      title('Water saturation - implicit solver')
      view([0, 0]), axis tight equal
      caxis([0 1]), colorbar      
           
      delete(h3);
      subplot(3,1,3),
      h3 = plotCellData(G, resSol_impl.s(:,1));
      title('Water saturation - theta scheme')
      view([0, 0]), axis tight equal
      caxis([0 1]), colorbar

   end

   % Update solution of pressure equation.
   resSol      = solveIncompFlow(resSol     , [], G, S, ...
                                 fluid, 'src', src);
   resSol_impl = solveIncompFlow(resSol_impl, [], G, S, ...
                                 fluid, 'src', src);
   resSol_t    = solveIncompFlow(resSol_t, [], G, S, ...
                                 fluid, 'src', src);

   drawnow;
   t = t + timestep;
end

normS(1) = norm(resSol.s(:,1) - resSol_t.s(:,1))/norm(resSol.s(:,1));
normS(2) = norm(resSol_impl.s(:,1) - resSol_t.s(:,1))/norm(resSol_impl.s(:,1));


