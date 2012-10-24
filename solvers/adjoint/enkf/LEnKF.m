
%----------------------------------------------------------------------------------
% SYNOPSIS:
%
% DESCRIPTION:
%   EnKF update with covariance localization in (x,y) space.
%
% PARAMETERS:
%
% RETURNS:
%   Updated ensemble of model states U
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
nx = grid.cartDims(1);
ny = grid.cartDims(2);
nz = grid.cartDims(3);
ne = nmembers;

U = A;

for ix = 1 : nx
for iy = 1 : ny
    
  if (iy == 1), disp(['ix = ' num2str(ix) ' (' num2str(nx) ')']), end

  % Find local well data
  lobs = []; dobs = [];
  id = find(loc > 0);
  if ~isempty(id)
      for j = 1 : length(id)
         [xj, yj, zj] = getIndex(loc(id(j)),nx,ny);
         do = sqrt((xj-ix)^2 + (yj-iy)^2);
         if do < range
           lobs = [lobs id(j)];
           dobs = [dobs do];
         end
      end
  end
  
  % Find local grid data
  id = find(loc < 0);
  if ~isempty(id)
      for i = max(1,ix - 2*range + 1) : min(nx,ix + 2*range - 1)
        for j = max(1,iy - 2*range + 1) : min(ny,iy + 2*range - 1)
          do = sqrt((i-ix)^2 + (j-iy)^2);
          if do < 2*range
            jd = ((1:nz)-1)*nx*ny + (j-1)*nx + i;
            [c, ia, iloc] = intersect(jd, abs(loc(id)));
            if ~isempty(iloc)
              lobs = [lobs iloc];
              dobs = [dobs do*ones(numel(iloc),1)];
            end
          end
        end
      end
  end

  if ~isempty(lobs)
  
  % find local state vector (currently only tested on 2D grids without non-active cells)
  Al = zeros(length(nstat),ne);
  id = (iy-1)*nx + ix;
  for i = 1 : length(nstat)
    var = A(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:);   
    Al(i,:) = var(id,:);
  end
  Am = mean(Al,2);
  Al = Al - repmat(Am,1,ne);    

%   Al = zeros(length(nstat)*nz,ne);
%   id = ((1:nz)-1)*nx*ny + (iy-1)*nx + ix;
%   for i = 1 : length(nstat)
%     var = A(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:);   
%     Al((i-1)*nz+1:i*nz,:) = var(id,:);
%   end
%   Am = mean(Al,2);
%   Al = Al - repmat(Am,1,ne);        
      
  % construct local matrices
  HA = D(lobs,:);
  Sl = (HA-repmat(mean(HA,2),1,ne))/sqrt(ne-1);
  Yl = Y(lobs,:);
  Rl = R(lobs);
  
  % update incorporating covariance localization
  rho = zeros(size(Al,1),length(lobs));
  for i = 1 : length(lobs)
      rho(:,i) = gaspari(dobs(i),range);
  end
  P1 = rho .* (Al * Sl');
  
  P2 = Sl*Sl';
  rho = eye(length(lobs));
  for i = 1 : length(lobs)
    [xi, yi, zi] = getIndex(abs(loc(lobs(i))),nx,ny);
    for j = i + 1 : length(lobs)
      [xj, yj, zj] = getIndex(abs(loc(lobs(j))),nx,ny);
      z = sqrt((xj-xi)^2 + (yj-yi)^2);
      rho(i,j) = gaspari(z,range);
      rho(j,i) = rho(i,j);
    end
  end

  P2 = rho .* P2 + diag(Rl);   
  W = P1 / P2;

  Ul = Al + repmat(Am,1,ne) + beta * sqrt(ne-1) .\ (W * (Yl - HA));

  for i = 1 : length(nstat)
    var = U(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:);
    var(id,:) = Ul((i-1)*nz+1:i*nz,:);
    U(sum(nstat(1:i-1))+1:sum(nstat(1:i)),:) = var;
  end

  end
  
end
end

clear lobs dobs Al Am HA Sl Yl Rl P1 P2 Ul W rho z
