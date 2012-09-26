function [s, varargout] = twophaseUpwReorder(s, tf, q, gm, pv, varargin)
%Single point upwind solver for Buckley-Leverett flow based on reordering.
%
% SYNOPSIS:
%   s = twophaseUpwReorder(s, t, q, gm, porvol, varargin)

%{
#COPYRIGHT#
%}

% $Id: twophaseUpwReorder.m 1953 2009-03-31 10:54:12Z bska $

opt = struct('mobility_ratio', 10.0,       ...
             'max_iterations', int32(100), ...  % nonlinear iterations
             'substeps',       int32(1),   ...
             'tolerance',      1e-6,       ...  % nonlinear tolerance
             'verbosity',      int32(2),   ...  % not used
             'scalar_solver',  'ridder');       % not used

opt = merge_options(opt, varargin{:});

opt.verbosity      = int32(opt.verbosity);
opt.max_iterations = int32(opt.max_iterations);
opt.substeps       = int32(opt.substeps);



if nargout == 1
    s       = implicitupwind(-gm', pv, s, tf, full(q), opt);
elseif nargout == 2
    [s, report] = implicitupwind(-gm', pv, s, tf, full(q), opt);
    varargout{1} = report;
else
    error('Too many output arguments');
end
