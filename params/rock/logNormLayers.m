function [K, L] = logNormLayers(N, M, varargin)
%Compute realization of lognormal, isotropic permeability field.
%
% SYNOPSIS:
%   K     = logNormLayers(N)
%   [K,L] = logNormLayers(N, M, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function creates a nx-by-ny-by-nz scalar permeability field
%   consisting of M layers. If M is a vector, then the number of layers
%   equals numel(M) and M gives the desired mean for each layer
%
% PARAMETERS:
%   N       - Three-element vector, [nx, ny, nz], specifying the number of
%             discrete values in the 'x', 'y', and 'z' coordinate
%             directions respectively.
%
%   M      -  The number of layers if a scalar and the desired means of the
%             layers if a vector.
%
%   'pn'/pv - List of 'key'/value pairs designating optional parameters.
%             Currently supported parameters are
%               - Verbose -- Whether or not to emit progress reports
%                            during the computation.
%                            Logical.  Default value: FALSE.
%
% RETURNS:
%   K - The scalar nx-by-ny-by-nz permeability field
%   L - run-length encoded list of layers in K

%{
#COPYRIGHT#
%}

% $Id: logNormLayers.m 1951 2009-03-31 10:28:32Z bska $

opt = struct('Verbose', false, ...
             'sigma',   3);
opt = merge_options(opt, varargin{:});

rand_layers = @(m,n) unique([1, fix(1 + rand([1,m-1])*(n-1)), n+1]);
% Process input parameters
if (nargin < 2) || isempty(M),
   L = [1:N(3), N(3)+1];  lmean = [];
elseif numel(M) == 1,
   % i.e., only the number of layers is specified
   if M > N(3),
      error(msgid('NumLayers:Excessive'), ...
           ['Number of requested permeability layers (%d) ', ...
            'exceeds number of model layers (%d).'], M, N(3));
   end
   if M < 0,
      L = [1, N(3)+1];
   elseif M == N(3),
      L = 1 : N(3) + 1;
   else
      L = [];
      while numel(L) < M + 1,
         L = rand_layers(M, N(3));
      end
   end
   lmean = [];
else
   % i.e., a mean is specified per layer
   if numel(M) > N(3),
      error(msgid('NumLayers:Excessive'), ...
           ['Number of requested permeability layers (%d) ', ...
            'exceeds number of model layers (%d).'], numel(M), N(3));
   end
   if numel(M) == N(3),
      L = 1 : N(3) + 1;
   else
      L = [];
      while numel(L) < numel(M) + 1,
         L = rand_layers(numel(M), N(3));
      end
   end
   lmean = M(:);
end

if opt.Verbose,
   fprintf('\nGenerating lognormal, layered permeability\n  Layers: [ ');
end
for i = 1 : length(L)-1,
    n = L(i+1) - L(i);
    k = smooth3(randn([N(1:2), n+2]) - 0.6*randn(1), ...
                'gaussian', [9, 3, 3], 4.5);
    k = exp(2 + opt.sigma*k);
    k = k(:,:,2:n+1);
    if numel(lmean) > 0,
       K(:,:,L(i):L(i+1)-1) = lmean(i)*k / mean(k(:));
    else
       K(:,:,L(i):L(i+1)-1) = k;
    end
    if opt.Verbose, fprintf('%d:%d ', L(i), L(i+1)-1); end
end
K = K(:);

if opt.Verbose,
   fprintf(']\n');
   fprintf('  min: %g, max: %g [mD], ratio: %g\n\n', ...
           min(K), max(K), max(K) / min(K));
end
