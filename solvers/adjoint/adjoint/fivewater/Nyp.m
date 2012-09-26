%
%     function [ mat ] = Nyp( y )
%
%     Purpose:
%
%     Given y compute the Jacobian N'(y) of the nonlinear term N(y) 
%     in Burgers equation
%
%     Parameters
%
%     On entry:
%
%     y      state y at one time step
%
%     On return:
%
%     mat    the matrix N'( y )  in sparse matrix format
%
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

  function [ mat ] = Nyp( y )
  
  ny = length(y);
  
  % subdiagonal
  e1         = zeros(ny,1);
  e1(1:ny-1) = -2*y(1:ny-1)-y(2:ny);
  
  % diagonal
  e2         = zeros(ny,1);
  e2(1)      = y(2);
  e2(2:ny-1) = y(3:ny)-y(1:ny-2);
  e2(ny)     = -y(ny-1);
  
  % superdiagonal
  e3         = zeros(ny,1);
  e3(2:ny)   = y(1:ny-1) + 2*y(2:ny);
  
  mat = (1/6)*spdiags([e1 e2 e3], -1:1, ny, ny);
      
%
% End of Npy.
