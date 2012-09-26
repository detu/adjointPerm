function M = computeElemMat3D()
%Calculate inner products of 3D basis functions for Taylor-Hood elements.
%
% SYNOPSIS:
%   M = computeElemMat3D()
%
% DESCRIPTION:
%   The velocity basis functions (bf) for each cell are numbered in the
%   following way:
%
%   7--8--9               16-17-18               25-26-27
%   |  |  | for the       |  |  |  for the       |  |  |  for the
%   4--5--6 lowest layer, 13-14-15 middle layer, 22-23-24 top layer.
%   |  |  |               |  |  |                |  |  |
%   1--2--3               10-11-12               19-20-21
%
%   The pressure basis functions for each cell are numbered in the
%   following way:
%
%   3-----4               7-----8
%   |  |  | for the       |  |  |  for the
%   |--|--|lowest layer,  |--|--| top layer.
%   |  |  |               |  |  |
%   1-----2               5-----6
%
% RETURNS:
%   M - System structure having the following fields:
%
%        - VV   - The inner product of the velocity basis functions (bf).
%        - VxVx - The inner product of the velocity bf differentiated with
%                 respect to x.
%        - VyVy - The inner product of the velocity bf differentiated with
%                 respect to y.
%        - VzVz - The inner product of the velocity bf differentiated with
%                 respect to z.
%        - VxP  - The inner product of the velocity bf differentiated with
%                 respect to x with the pressure bf.
%        - VyP  - The inner product of the velocity bf differentiated with
%                 respect to y with the pressure bf.
%        - VzP  - The inner product of the velocity bf differentiated with
%                 respect to z with the pressure bf.
%
% COMMENTS:
%   The matrices are stored permanently in the files ElementMat3D.mat and
%   BasisFunc3D.mat in the directory
%
%      fullfile(ROOTDIR, 'solvers', 'stokes-brinkman')

%{
#COPYRIGHT#
%}

% $Date: 2009-10-01 17:08:49 +0200 (to, 01 okt 2009) $
% $Revision: 2930 $

% Read basis functions
   [W, dWdx, dWdy, dWdz, L] = basisfunc3D();

   [nVdof1D, nPdof1D, dim] = deal(3, 2, 3);

   % The inner product of the velocity basis functions and their derivates.
   [VV, VxVx, VyVy, VzVz] = deal(zeros([nVdof1D, nVdof1D] .^ dim));
   for i = 1 : nVdof1D ^ dim,
      Vi  = W{i};
      Vxi = dWdx{i};
      Vyi = dWdy{i};
      Vzi = dWdz{i};

      for j = 1 : nVdof1D ^ dim,
         Vj  = W{j};
         Vxj = dWdx{j};
         Vyj = dWdy{j};
         Vzj = dWdz{j};

         integrand1 = @(x,y,z) Vi (x,y,z) .* Vj (x,y,z);
         VV(i,j)    = triplequad(integrand1, -1, 1, -1, 1, -1, 1);

         integrand2 = @(x,y,z) Vxi(x,y,z) .* Vxj(x,y,z);
         VxVx(i,j)  = triplequad(integrand2, -1, 1, -1, 1, -1, 1);

         integrand3 = @(x,y,z) Vyi(x,y,z) .* Vyj(x,y,z);
         VyVy(i,j)  = triplequad(integrand3, -1, 1, -1, 1, -1, 1);

         integrand4 = @(x,y,z) Vzi(x,y,z) .* Vzj(x,y,z);
         VzVz(i,j)  = triplequad(integrand4, -1, 1, -1, 1, -1, 1);
      end
   end

   % The inner product of the derivates of the velocity basis functions
   % with the pressure basis functions.
   [VxP, VyP, VzP] = deal(zeros([nVdof1D, nPdof1D] .^ dim));
   for i = 1 : nVdof1D ^ dim,
      Vxi = dWdx{i};
      Vyi = dWdy{i};
      Vzi = dWdz{i};

      for j = 1 : nPdof1D ^ dim,
         Lj = L{j};

         integrand1 = @(x,y,z) Vxi(x,y,z) .* Lj(x,y,z);
         VxP(i,j)   = triplequad(integrand1, -1, 1, -1, 1, -1, 1);

         integrand2 = @(x,y,z) Vyi(x,y,z) .* Lj(x,y,z);
         VyP(i,j)   = triplequad(integrand2, -1, 1, -1, 1, -1, 1);

         integrand3 = @(x,y,z) Vzi(x,y,z) .* Lj(x,y,z);
         VzP(i,j)   = triplequad(integrand3, -1, 1, -1, 1, -1, 1);
      end
   end

   M.VV   = VV  ;
   M.VxVx = VxVx;   M.VyVy = VyVy;   M.VzVz = VzVz;
   M.VxP  = VxP ;   M.VyP  = VyP ;   M.VzP  = VzP ;

   save(sbpath('ElementMat3D'), 'M');
end

%--------------------------------------------------------------------------

function [W, dWdx, dWdy, dWdz, L] = basisfunc3D()
% basisfunc - Calculates the basis functions
%
% SYNOPSIS:
%   [W, dWdx, dWdy, dWdz, L] = basisfunc3D()
%
% RETURNS:
%   W    - The velocity basis functions.
%
%   dWdx - Velocity basis functions differentiated with respect to x.
%
%   dWdy - Velocity basis functions differentiated with respect to y.
%
%   dWdz - Velocity basis functions differentiated with respect to z.
%
%   L    - The pressure basis functions.

   [nVdof1D, nPdof1D, dim] = deal(3, 2, 3);

   % The velocity basis functions
   B1 = @(x) -0.5.*x.*(1-x);
   B2 = @(x)  (1-x.^2);
   B3 = @(x)  0.5.*x.*(1+x);
   B  = {B1, B2, B3};

   W  = cell([1, nVdof1D ^ dim]);
   ix = 0;
   for k = 1 : nVdof1D,       BK = B{k};
      for j = 1 : nVdof1D,    BJ = B{j};
         for i = 1 : nVdof1D, BI = B{i};
            ix = ix + 1;
            W{ix} = @(x,y,z) BI(x) .* BJ(y) .* BK(z);
         end
      end
   end

   % The velocity basis functions derived by x
   dB1 = @(x) -0.5*(1-2*x);
   dB2 = @(x) -2*x;
   dB3 = @(x)  0.5*(1+2*x);
   dB  = {dB1, dB2, dB3};

   dWdx = cell([1, nVdof1D ^ dim]);
   ix   = 0;
   for k = 1 : nVdof1D,       BK =  B{k};
      for j = 1 : nVdof1D,    BJ =  B{j};
         for i = 1 : nVdof1D, BI = dB{i};
            ix = ix + 1;
            dWdx{ix} = @(x,y,z) BI(x) .* BJ(y) .* BK(z);
         end
      end
   end

   % The velocity basis functions derived by y
   dWdy = cell([1, nVdof1D ^ dim]);
   ix   = 0;
   for k = 1 : nVdof1D,       BK =  B{k};
      for j = 1 : nVdof1D,    BJ = dB{j};
         for i = 1 : nVdof1D, BI =  B{i};
            ix = ix + 1;
            dWdy{ix} = @(x,y,z) BI(x) .* BJ(y) .* BK(z);
         end
      end
   end

   % The velocity basis functions derived by z
   dWdz = cell([1, nVdof1D ^ dim]);
   ix = 0;
   for k = 1 : nVdof1D,       BK = dB{k};
      for j = 1 : nVdof1D,    BJ =  B{j};
         for i = 1 : nVdof1D, BI =  B{i};
            ix = ix + 1;
            dWdz{ix} = @(x,y,z) BI(x) .* BJ(y) .* BK(z);
         end
      end
   end

   % The pressure basis functions
   B1 = @(x) 0.5*(1-x);
   B2 = @(x) 0.5*(1+x);
   Bl = {B1, B2};

   L = cell([1, nPdof1D ^ dim]);
   ix = 0;
   for k = 1 : nPdof1D,       BK = Bl{k};
      for j = 1 : nPdof1D,    BJ = Bl{j};
         for i = 1 : nPdof1D, BI = Bl{i};
            ix = ix + 1;
            L{ix} = @(x,y,z) BI(x) .* BJ(y) .* BK(z);
         end
      end
   end

   % The integrated velocity basis functions
   [W_1yz, W1yz, Wx_1z, Wx1z, Wxy_1, Wxy1] = deal(cell([1, nVdof1D ^ dim]));
   for j = 1 : nVdof1D ^ dim,
      W_1yz{j} = dblquad(@(  y,z) W{j}(-1, y, z), -1, 1, -1, 1);
      W1yz {j} = dblquad(@(  y,z) W{j}( 1, y, z), -1, 1, -1, 1);
      Wx_1z{j} = dblquad(@(x,  z) W{j}( x,-1, z), -1, 1, -1, 1);
      Wx1z {j} = dblquad(@(x,  z) W{j}( x, 1, z), -1, 1, -1, 1);
      Wxy_1{j} = dblquad(@(x,y  ) W{j}( x, y,-1), -1, 1, -1, 1);
      Wxy1 {j} = dblquad(@(x,y  ) W{j}( x, y, 1), -1, 1, -1, 1);
   end

   save(sbpath('BasisFunc3D'), 'W', 'dWdx', 'dWdy', 'dWdz', 'L', ...
        'W_1yz', 'W1yz', 'Wx_1z', 'Wx1z', 'Wxy_1', 'Wxy1');
end

%--------------------------------------------------------------------------

function p = sbpath(p)
   p = fullfile(ROOTDIR, 'solvers', 'stokes-brinkman', p);
end
