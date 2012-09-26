function sim = readEclipseFile(filename)
%Input ECLIPSE or FrontSim Deck specification

%  param = readPARAMS(filename);
%  param = append_fields(param, readGRDECL(filename));
param = readGRDECL(filename);

if ~all(isfield(param, {'cartDims'}))
  error('Error in input file');
end

sim.grid    = getGrid(param);
sim.rock    = getRock(param);
sim.ressol  = getInitialData(param);
sim.bc      = getBoundaryConditions(param);
sim.relperm = getRelPerm(param);
sim.pvt     = getPVT(param);
sim.dt      = getfield(param, 'TSTEP');
sim.param   = param;

clf
plotGrid(sim.grid);view(30,30)


%%% =======================================================================
%   Private functions follow below
%

%%% -----------------------------------------------------------------------
function param = append_fields(param, fields) %#ok
fn = fieldnames(fields);
for n = 1 : numel(fn),
  param.(fn{n}) = fields.(fn{n});
end


%%% -----------------------------------------------------------------------
function bc = getBoundaryConditions(param)
% Boundary conditions are per default no-flow.
bc = struct;
if isfield(param, 'PSIDE')
  % Add Dirichlet boundary condition
end

if isfield(param, 'PSIDEH')
  % Add Dirichlet boundary condition
end

if isfield(param, 'FLUXSIDE')
  % Add Neumann boundary conditions
end

%%% -----------------------------------------------------------------------
function resSol = getInitialData(param)
resSol = struct;
n = prod(param.cartDims);

if isfield(param, 'SWAT')
  % Read inital water saturation
else
  resSol.sw = zeros(n,1);
end

if isfield(param, 'SGAS')
    % Read initial gas saturation
else
 resSol.sg = zeros(n,1);
end


if isfield(param, 'RS')
  % Read inital dissolved gas
else
  resSol.rl = zeros(n,1);
end

if isfield(param, 'RV')
  % Read inital vaporized oil
else
  resSol.rv = zeros(n,1);
end

%{
  if isfield(param, 'PRESSURE')
    % Read inital pressure
  elseif all(isfield(param, {'EQUIL', 'DENSITY'}))
    %
  else
    error ('Unable to set initial pressure');
  end

%}


%%% -----------------------------------------------------------------------
function relperm = getRelPerm(param) %#ok
relperm = struct;



%%% -----------------------------------------------------------------------
function pvt = getPVT(param)  %#ok
pvt = struct;



%%% -----------------------------------------------------------------------
function G = getGrid(param)
%
%  three types of grids are currently supported
%    - Tensor grids, vertical pillars and no faults as
%      specified by DXV,DXV,DXZ[,DEPTHZ] grids
%    - Matching hexahedral grids as specified by
%      COORD[XYZ]. NB! Cells may degenerate
%    - General Corner-point grids
G = struct;

if all(isfield(param, {'DXV', 'DYV', 'DZV', 'DEPTHZ'}))
  % Tensor grid: vertical pillars, no faults.
  if isfield(param,'ACTNUM')
    warning(100, 'ACTNUM not supported for tensor grids yet');
  end
  G = tensorGrid(param.DXV, param.DYV, param.DZV, param.DEPTHZ);

elseif all(isfield(param, {'COORDX', 'COORDY', 'COORDZ'}))
  % General matching grid given by node coordinates.
  G = buildCoordGrid(param);
  % Do something

elseif all(isfield(param, {'COORD', 'ZCORN'}))
  % Pillar grid
  G = buildMatchingGrid(param);
  % Do something

end

G = computeGeometry(G);


%%% -----------------------------------------------------------------------
function rock = getRock(param)
%
% rock.perm is [nx3](diagonal) or [nx6](full, symmetric)
% rock.poro is [n]
n = prod(param.cartDims);

if all(~isfield(param,{'PERMXY','PERMXZ','PERMYX',...
                      'PERMYZ','PERMZX','PERMZY'}))
  % diagonal permeability tensor
  if isfield(param, {'PERMX', 'PERMY', 'PERMZ'})
    rock.perm = zeros(n,3);
    rock.perm(:,1) = param.PERMX;
    rock.perm(:,2) = param.PERMY;
    rock.perm(:,3) = param.PERMZ;
  else
    error(['The permeability must at least be diagonal!'...
           'Specify PERMX, PERMY and PERMZ in parameter file.']);
  end
else
  %% full symmetric permeability tensor
  rock.perm = zeros(n, 6);
  if isfield(param, {'PERMX', 'PERMY', 'PERMZ'})
    rock.perm(:,1)=param.PERMX;
    rock.perm(:,4)=param.PERMY;
    rock.perm(:,6)=param.PERMZ;
  else
    error(['The permeability must at least be diagonal!'...
           'Specify PERMX, PERMY and PERMZ in parameter file.']);
  end

  if all(isfield(param, {'PERMXY','PERMYX'})) && ...
          any(param.PERMXY~=param.PERMYX),
    warning(['PERMXY and PERMYX differ.' ...
             'Asymmetric tensors are currently not supported.']); %#ok
  end
  if all(isfield(param, {'PERMXZ','PERMZX'})) && ...
      any(param.PERMXZ~=param.PERMZX),
    warning(['PERMXZ and PERMZX differ.' ...
             'Asymmetric tensors are currently not supported.']); %#ok
  end
  if all(isfield(param, {'PERMYZ','PERMZY'})) && ...
          any(param.PERMYZ~=param.PERMZY),
    warning(['PERMYZ and PERMZY differ.' ...
             'Asymmetric tensors are currently not supported.']); %#ok
  end

  if isfield(param, 'PERMXY'), rock.perm(:,2)=param.PERMXY; end
  if isfield(param, 'PERMXZ'), rock.perm(:,3)=param.PERMXZ; end
  if isfield(param, 'PERMYZ'), rock.perm(:,5)=param.PERMYZ; end
  if isfield(param, 'PERMYX'), rock.perm(:,2)=param.PERMYX; end
  if isfield(param, 'PERMZX'), rock.perm(:,3)=param.PERMZX; end
  if isfield(param, 'PERMZY'), rock.perm(:,5)=param.PERMZY; end
end

if isfield(param, 'PORO')
  rock.poro = param.PORO;
else
  error ('Porosity is missing!');
end
