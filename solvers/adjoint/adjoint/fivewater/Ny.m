%
%     function [ val ] = Ny( y )
%
%     Purpose:
%
%     Given y evaluate the nonlinear term N(y) in Burgers equation
%
%     Parameters
%
%     On entry:
%
%     y      state y at one time step
%
%     On return:
%
%     val    N( y )  
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%

  function [ val ] = Ny( y )
  
  ny = length(y);
  val = zeros(size(y)); 
  
  val(1) = y(1)*y(2) + y(2)^2;
  for i = 2:ny-1
     val(i) = -y(i-1)^2 - y(i)*y(i-1) + y(i)*y(i+1) + y(i+1)^2;
  end
  val(ny) = -y(ny-1)^2 - y(ny)*y(ny-1);
  
  val = val/6;
  
%
% End of Ny.
