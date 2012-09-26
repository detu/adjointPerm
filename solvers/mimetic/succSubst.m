function [xr, xw] = succSubst(xr, xw, solver, varargin)
%Successive substitution algorithm for compressible pressure system.
%
% SYNOPSIS:
%   [resSol, wellSol] = succSubst(resSol, wellSol, solver)
%   [resSol, wellSol] = succSubst(resSol, wellSol, solver, 'pn1', pv1, ...)
%
% PARAMETERS:
%   resSol  - Reservoir solution structure representing reservoir state at
%             previous time step.
%
%   wellSol - Well solution structure representing well state at previous
%             time step.
%
%   solver  - Pressure solver callback which implements the body of single
%             successive substitution algorithm.  Assumed to support the
%             calling syntax
%
%                  [resSol, wellSol] = solver(resSol, wellSol)
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Tol     -- Solver tolerance.  The process is terminated
%                            when NORM(dv,INF) < Tol with 'dv' denoting the
%                            velocity/flux increment.
%                            Double.  Default value: Tol = 5.0e-9.
%
%               - MaxIt   -- Maximum number of successive substitution
%                            iterations.  The process is terminated when
%                            'MaxIt' iterations have been executed.
%                            Integer.  Default value: MaxIt = 10.
%
%               - Verbose -- Whether or not to emit informational messages
%                            concerning the convergence history of the
%                            process.
%                            Logical.  Default value: Verbose = False.
%
% RETURNS:
%   resSol  - Updated reservoir solution structure.
%   wellSol - Updated well solution structure.
%
% EXAMPLE:
%   resSol  = initResSol (G, initP);
%   wellSol = initWellSol(W, initP);
%   tab     = readpvt(fullfile(SAMSIMROOTDIR, 'fluid', ...
%                              'Data', 'compressible.txt'));
%   DT = 1;
%   p0 = resSol.cellPressure;
%   solver = @(xr,xw) solveBlackOilWellSystem(xr, xw, G, rock, S, W, ...
%                                             tab, p0, DT, 'bc', bc, ...
%                                             'src', src);
%   [resSol, wellSol] = succSubst(resSol, wellSol, solver, ...
%                                 'MaxIt', 40, 'Verbose', true);
%   plotCellData(G, resSol.cellPressure)
%
% SEE ALSO:
%   solveBlackOilWellSystem, solveCoarsePSysBO.

%{
#COPYRIGHT#
%}

% $Id: succSubst.m 1953 2009-03-31 10:54:12Z bska $

   opt = struct('Tol', 5.0e-9, 'MaxIt', 10, 'Verbose', false);
   opt = merge_options(opt, varargin{:});

   opt.Tol   = abs(opt.Tol);    assert (opt.Tol   >  0);
   opt.MaxIt = abs(opt.MaxIt);  assert (opt.MaxIt >= 0);

   if opt.Verbose,
      fprintf('Solving pressure by successive substitution...\n')
      nDigits = floor(log10(opt.MaxIt)) + 1;
      tic
   end

   it = 0;  res = 10 * max(opt.Tol, 1);
   while (res > opt.Tol) && (it < opt.MaxIt),
      v0 = xr.cellFlux;

      [xr, xw] = solver(xr, xw);

      res = norm(xr.cellFlux - v0, inf);

      it = it + 1;

      if opt.Verbose,
         fprintf('   %0*d: Norm(dv,inf) = %10.4e\n', nDigits, it, res)
      end
   end

   if opt.Verbose, toc, end
   if res > opt.Tol,
      warning(msgid('NlTol:NotMet'), ...
              '   Non-linear criterion not met: %10.4e > %10.4e\n', ...
              res, opt.Tol);
   end
end
