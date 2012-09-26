function [J, grad] = coarsegrad()
% This is an example Matlab code that reads Eclipse summary and restart
% files and calculated objective function and gradients.
%
% Sigurd Ivar Aanonsen, CIPR
% March 2008

%********************************************
% Load simulated and measured data
%********************************************

simhistfile     = sprintf('ECLIPSE.UNSMRY');
filecontents_sh = readXfile(simhistfile,0);
simhist         = filecontents_sh.PARAMS;
histfile        = sprintf('REFCASE.UNSMRY');
filecontents_h  = readXfile(histfile,0);
hist            = filecontents_h.PARAMS;
[n_dsim,n_tsim] = size(simhist);
[n_d,n_t]       = size(hist);
if n_d ~= n_dsim
    fprintf(1, 'n_d    = %d\n', n_d);
    fprintf(1, 'n_dsim = %d\n', n_dsim);
    error ('Not the same no of data in hist & sim files')
end
if n_t ~= n_tsim
    fprintf(1, 'n_t    = %d\n', n_t);
    fprintf(1, 'n_tsim = %d\n', n_tsim);
    error ('Not the same no of time steps in hist & sim files')
end
%'No of data:'
n_d = n_d - 2;
fprintf(1, 'n_d    = %d\n', n_d);
%'No of time steps:'
fprintf(1, 'n_t    = %d\n', n_t);

% load start.mat
%********************************************
% Calculate measure error as one percent
% of the pressure difference between the
% production well pressure and the maximum
% injection well pressure in reference data
%********************************************
pbhp       = 2000.0; % Production well pressure
vmax       = max(hist);
pmax       = max(vmax);
pdiff      = pmax - pbhp;
pn         = 1;
meas_error = pdiff*pn/100;
%********************************************
% Add 1% Gaussian random noise to data
%********************************************
for i = 1:n_d
    for j = 1:n_t
        ibhp           = simhist(i+2,j);
        simhist(i+2,j) = ibhp + meas_error*randn;
    end
end

%********************************************
% Calculate and save objective function
%********************************************

objfun = 0;
r      = zeros(n_d, n_t);
for i = 1:n_d
  for j = 1:n_t
      r(i,j) = simhist(i+2,j) - hist(i+2,j);
      objfun = objfun + ( r(i,j)/meas_error )^2;
  end
end
rms   = n_d*n_t-npp + sqrt(2*(n_d*n_t-npp));
xx(1) = double(objfun);
xx(2) = double(rms);
save objfunct.dat xx -ASCII;

%********************************************
% Load gradients
%********************************************

%'Reading gradients:'
for j = 1:n_t
    if j < 10
        xfile = strcat('ECLIPSE.X000',sprintf('%d',j));
    elseif j < 100
        xfile = strcat('ECLIPSE.X00',sprintf('%d',j));
    elseif j < 1000
        xfile = strcat('ECLIPSE.X0',sprintf('%d',j));
    else
        error('Too many time steps')
    end

    filecontents  = readXfile(xfile,0);
    ajgpara       = filecontents.AJGPARA;
    if j == 1
        n_p = ajgpara(4);
        fprintf(1, 'n_p    = %d\n', n_p);
    end
    if n_p ~= ajgpara(4);
        fprintf(1, 'n_p    = %d\n', n_p);
        fprintf(1, 'n_d    = %f\n', ajgpara(4));
        error ('Wrong no of parameters in gradients file')
    end
    ajgnames  = filecontents.AJGNAMES;
    n_g = length(ajgnames)/4;

    try
      i = 1;
      temp=filecontents.AJGFN1;
      ajgfn(:,i,j)=temp;
      i = 2;
      temp=filecontents.AJGFN2;
      ajgfn(:,i,j)=temp;
      i = 3;
      temp=filecontents.AJGFN3;
      ajgfn(:,i,j)=temp;
      i = 4;
      temp=filecontents.AJGFN4;
      ajgfn(:,i,j)=temp;
      i = 5;
      temp=filecontents.AJGFN5;
      ajgfn(:,i,j)=temp;
      i = 6;
      temp=filecontents.AJGFN6;
      ajgfn(:,i,j)=temp;
      i = 7;
      temp=filecontents.AJGFN7;
      ajgfn(:,i,j)=temp;
      i = 8;
      temp=filecontents.AJGFN8;
      ajgfn(:,i,j)=temp;
      i = 9;
      temp=filecontents.AJGFN9;
      ajgfn(:,i,j)=temp;
      i = 10;
      temp=filecontents.AJGFN10;
      ajgfn(:,i,j)=temp;
      i = 11;
      temp=filecontents.AJGFN11;
      ajgfn(:,i,j)=temp;
      i = 12;
      temp=filecontents.AJGFN12;
      ajgfn(:,i,j)=temp;
      i = 13;
      temp=filecontents.AJGFN13;
      ajgfn(:,i,j)=temp;
      i = 14;
      temp=filecontents.AJGFN14;
      ajgfn(:,i,j)=temp;
      i = 15;
      temp=filecontents.AJGFN15;
      ajgfn(:,i,j)=temp;
      i = 16;
      temp=filecontents.AJGFN16;
      ajgfn(:,i,j)=temp;
      i = 17;
      temp=filecontents.AJGFN17;
      error ('Not more than 16 datatypes allowed')
    catch exception
      rethrow(exception);  
    end
end

%********************************************
% Calculate and save gradient of objective function
%********************************************
objgrad = zeros(1,n_p);
% Fine scale gradient
for k = 1:n_p  % Parameter 1 (PERMXY)
    for j = 1:n_t
        for i = 1:n_d
            grad       = ajgfn(k,i,j);
            objgrad(k) = objgrad(k) + r(i,j)*grad;
        end
        objgrad(k) = objgrad(k)/(meas_error^2);
    end
end

grad = objgrad;
J    = objfun;
