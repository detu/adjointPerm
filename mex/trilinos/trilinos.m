function x = trilinos(A, b)
% trilinos -- Wrap multi-level linear solver library Trilinos.
%
% SYNOPSIS:
%   x = trilinos(A, b)
%
% DESCRIPTION:
%   Solves the system Ax = b of linear equations by means of the Trilinos
%   multi-level solver package distributed by Sandia National Labs:
%
%       http://trilinos.sandia.gov/
%
% PARAMETERS:
%   A - System matrix.  Assumed to be SPARSE.
%   b - System right hand side.
%
% RETURNS:
%   x - System solution.

% $Id: trilinos.m 1738 2009-03-17 09:15:50Z jrn $

if false
   fn = tempname('.');
   [f, msg] = fopen (fn, 'wt');
   if f < 0, error(msg); end
   assert(issparse(A));
   fprintf(f, ['/estimate_condition_number=false\n',    ...
               '/maximum_iterations=1000\n',            ...
               '/solver_library=trilinos\n',            ...
               '/solver_type=CG\n',                     ...
               '/tolerance=1e-10\n',                    ...
               '/use_multilevel_preconditioner=true\n', ...
               '/x_result_do_copy=true\n',              ...
               '/solver_output=10\n']);
   fclose(f);

   c = onCleanup(@() delete(fn)); %#ok
   
else
   fn = fullfile(ROOTDIR, 'mex','trilinos', 'trilinos_cg_ml.xml');
end
try
   x  = matlab_solve(A, b, fn);
catch ME
   idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
   if strcmp(idSegLast, 'UndefinedFunction')
      error('In trilinos.m:  mex function matlab_solve.mex* not avaliable');
   end
   rethrow(ME);
end
