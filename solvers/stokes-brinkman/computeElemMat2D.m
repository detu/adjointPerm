function M = computeElemMat2D()
%Calculate inner products of 2D basis functions for Taylor-Hood elements.
%
% SYNOPSIS:
%   M = computeElemMat2D()
%
% DESCRIPTION:
%   The velocity basis functions for each cell are numbered in the
%   following way:
%
%   7--8--9
%   |  |  |
%   4--5--6
%   |  |  |
%   1--2--3
%
%   The pressure basis functions for each cell are numbered in the
%   following way:
%
%   3-----4
%   |  |  |
%   |--|--|
%   |  |  |
%   1-----2
%
% RETURNS:
%   M - System structure having the following fields:
%
%        - VV   - The inner product of the velocity basis functions (bf).
%        - VxVx - The inner product of the velocity bf differentiated with
%                 respect to x.
%        - VyVy - The inner product of the velocity bf differentiated with
%                 respeect to y
%        - VxP  - The inner product of the velocity bf differentiated with
%                 respect to x with the pressure bf.
%        - VyP  - The inner product of the velocity bf differentiated with
%                 respect to y with the pressure bf.
%
% COMMENTS:
%   The matrices are stored permanently in the files ElementMat2D.mat and
%   BasisFunc2D.mat in the directory
%
%      fullfile(ROOTDIR, 'solvers', 'stokes-brinkman')

% $Date: 2009-10-01 17:08:49 +0200 (to, 01 okt 2009) $
% $Revision: 2930 $

   % Read basis functions
   [W, dWdx, dWdy, L] = basisfunc2D();

   [nVdof1D, nPdof1D, dim] = deal(3, 2, 2);

   % The inner product of the velocity basis functions and their derivates.
   [VV, VxVx, VyVy] = deal(zeros([nVdof1D, nPdof1D] .^ dim));
   for i = 1 : nVdof1D ^ dim,
      Vi  = W{i};
      Vxi = dWdx{i};
      Vyi = dWdy{i};

      for j = 1 : nVdof1D ^ dim,
         Vj  = W{j};
         Vxj = dWdx{j};
         Vyj = dWdy{j};

         integrand1 = @(x,y) Vi (x,y) .* Vj (x,y);
         VV(i,j)    = dblquad(integrand1, -1, 1, -1, 1);

         integrand2 = @(x,y) Vxi(x,y) .* Vxj(x,y);
         VxVx(i,j)  = dblquad(integrand2, -1, 1, -1, 1);

         integrand3 = @(x,y) Vyi(x,y) .* Vyj(x,y);
         VyVy(i,j)  = dblquad(integrand3, -1, 1, -1, 1);
      end
   end

   % The inner product of the derivates of the velocity basis functions
   % with the pressure basis functions.
   [VxP, VyP] = deal(zeros([nVdof1D, nPdof1D] .^ dim));
   for i = 1 : nVdof1D ^ dim,
      Vxi = dWdx{i};
      Vyi = dWdy{i};

      for j = 1 : nPdof1D ^ dim,
         Lj = L{j};

         integrand1 = @(x,y) Vxi(x,y) .* Lj(x,y);
         VxP(i,j)   = dblquad(integrand1, -1, 1, -1, 1);

         integrand2 = @(x,y) Vyi(x,y) .* Lj(x,y);
         VyP(i,j)   = dblquad(integrand2, -1, 1, -1, 1);
      end
   end

   % The inner product of the velocity basis functions and their derivates.
   VVx = zeros([nVdof1D, nVdof1D] .^ dim);
   for i = 1 : nVdof1D ^ dim,
      Vi = W{i};
      for j = 1 : nVdof1D ^ dim,
         Vxj = dWdx{j};

         integrand = @(x,y) Vi(x,y) .* Vxj(x,y);
         VVx(i,j)  = dblquad(integrand, -1, 1, -1, 1);
      end
   end

   M.VV   = VV  ;
   M.VxVx = VxVx;   M.VyVy = VyVy;
   M.VxP  = VxP ;   M.VyP  = VyP ;

   save(sbpath('ElementMat2D'), 'M');
end

%--------------------------------------------------------------------------

function [W, dWdx, dWdy, L] = basisfunc2D()
% basisfunc - Calculates the basis functions
%
% SYNOPSIS:
%   [W, dWdx, dWdy, L]=basisfunc2D()
%
% RETURNS:
%   W    - The velocity basis functions
%
%   dWdx - The velocity basis functions derived by x
%
%   dWdy - The velocity basis functions derived by y
%
%   L    - The pressure basis functions

   % Velocity basis functions.
   W    = cell([1, 9]);
   W{1} = @(x,y)  x.*y./4.*(1-x).*(1-y);
   W{2} = @(x,y) -y./2.*(1-x.^2).*(1-y);
   W{3} = @(x,y) -x.*y./4.*(1+x).*(1-y);
   W{4} = @(x,y) -x./2.*(1-x).*(1-y.^2);
   W{5} = @(x,y) (1-x.^2).*(1-y.^2);
   W{6} = @(x,y)  x./2.*(1+x).*(1-y.^2);
   W{7} = @(x,y) -x.*y./4.*(1-x).*(1+y);
   W{8} = @(x,y)  y./2.*(1-x.^2).*(1+y);
   W{9} = @(x,y)  x.*y./4.*(1+x).*(1+y);

   % Velocity basis functions differentiated with respect to x.
   dWdx    = cell([1, 9]);
   dWdx{1} = @(x,y)  y./4.*(1-y).*(1-2.*x);
   dWdx{2} = @(x,y)  x.*y.*(1-y);
   dWdx{3} = @(x,y) -1./4.*y.*(1+2.*x).*(1-y);
   dWdx{4} = @(x,y) -1./2.*(1-2.*x).*(1-y.^2);
   dWdx{5} = @(x,y) -2.*x.*(1-y.^2);
   dWdx{6} = @(x,y)  1./2.*(1+2.*x).*(1-y.^2);
   dWdx{7} = @(x,y) -1./4.*y.*(1-2.*x).*(1+y);
   dWdx{8} = @(x,y)  -x.*y.*(1+y);
   dWdx{9} = @(x,y)  1./4.*y.*(1+2.*x).*(1+y);

   % Velocity basis functions differentiated with respect to y.
   dWdy    = cell([1, 9]);
   dWdy{1} = @(x,y)  x./4.*(1-x).*(1-2.*y);
   dWdy{2} = @(x,y) -1./2.*(1-x.^2).*(1-2.*y);
   dWdy{3} = @(x,y) -1./4.*x.*(1+x).*(1-2.*y);
   dWdy{4} = @(x,y)  x.*(1-x).*y;
   dWdy{5} = @(x,y) -2.*y.*(1-x.^2);
   dWdy{6} = @(x,y) -x.*(1+x).*y;
   dWdy{7} = @(x,y) -1./4.*x.*(1-x).*(1+2.*y);
   dWdy{8} = @(x,y)  1./2.*(1-x.^2).*(1+2.*y);
   dWdy{9} = @(x,y)  1./4.*x.*(1+x).*(1+2.*y);

   % Pressure basis functions.
   L1=@(x,y) (1-x).*(1-y)./4;   % == L4(-x,-y)
   L2=@(x,y) (1+x).*(1-y)./4;   % == L4( x,-y)
   L3=@(x,y) (1-x).*(1+y)./4;   % == L4(-x, y)
   L4=@(x,y) (1+x).*(1+y)./4;

   L ={L1  L2  L3  L4};

   % Integrated velocity basis functions.
   [W_1y, W1y, Wx_1, Wx1] = deal(cell([1, 9]));
   for j = 1 : 9,
      W_1y{j} = quad(@(y) W{j}(-1, y), -1, 1);
      W1y {j} = quad(@(y) W{j}( 1, y), -1, 1);
      Wx_1{j} = quad(@(x) W{j}( x,-1), -1, 1);
      Wx1 {j} = quad(@(x) W{j}( x, 1), -1, 1);
   end

   save(sbpath('BasisFunc2D'), 'W', 'dWdx', 'dWdy', 'L', ...
        'W_1y', 'W1y', 'Wx_1', 'Wx1');
end

%--------------------------------------------------------------------------

function p = sbpath(p)
   p = fullfile(ROOTDIR, 'solvers', 'stokes-brinkman', p);
end
