function [f, df] = mismatchHM(x)

% mismatchHM.m This function returns the function value
% of mismatch reference model (measurement) with actual model, 
% and its first order derivative

% open PERMX.DAT file
% write it with 1 mD (initial values)
nx  = length(x);
fid = fopen('PERMX.DAT', 'w');
for k = 1:nx
    fprintf(fid, '%f\n', x(k));
end
fprintf(fid, '//');
fclose(fid);

% simulation call, i.e ECLIPSE300
fprintf(1,'************************************************\n');
fprintf(1,'Making a simulation run \n');
fprintf(1,'************************************************\n');
cmdline = ('/usr/local/ecl/macros/@e300 < ../code/ecl.in > ecl.out');
unix(cmdline);

% read restart and summary files and write required data to file
% compute objective function and gradient
fprintf(1,'************************************************');
fprintf(1,'Calling matlab to process Eclipse results');
fprintf(1,'************************************************');
[J, grad] = coarsegrad;
f         = J;

if nargout > 1
  df = grad;
end

