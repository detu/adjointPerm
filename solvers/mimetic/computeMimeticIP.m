function S = computeMimeticIP(G, rock, varargin)
%Compute mimetic inner product matrices.
%
% SYNOPSIS:
%   S = computeMimeticIP(G, rock)
%   S = computeMimeticIP(G, rock, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and FOUR/SIX columns for
%             a symmetric full tensor permeability.  In the latter case,
%             each cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Type -- The kind of system to assemble.  The choice made
%                         for this option influences which pressure solvers
%                         can be employed later.
%                         String.  Default value = 'hybrid'.
%                 Supported values are:
%                   - 'hybrid'      : Hybrid system with inverse of B
%                                     (for Schur complement reduction)
%                   - 'mixed'       : Hybrid system with B
%                                     (for reduction to mixed form)
%                   - 'tpfa'        : Hybrid system with B
%                                     (for reduction to tpfa form)
%                   - 'comp_hybrid' : Both 'hybrid' and 'mixed'
%
%               - InnerProduct -- The inner product with which to define
%                                 the mass matrix.
%                                 String.  Default value = 'ip_simple'.
%                 Supported values are:
%                   - 'ip_simple'   : inner product having the 2*tr(K)-term.
%                   - 'ip_tpf'      : inner product giving method equivalent
%                                     to two-point flux approximation (TPFA).
%                   - 'ip_quasitpf' : inner product ''close to'' TPFA
%                                     (equivalent for Cartesian grids with
%                                      diagonal permeability tensor).
%                   - 'ip_rt'       : Raviart-Thomas for Cartesian grids,
%                                     else not valid.
%
%               - Verbose -- Whether or not to emit progress reports during
%                            the assembly process.
%                            Logical.  Default value = FALSE.
%
% RETURNS:
%   S - Pressure linear system structure having the following fields:
%         - BI / B   : inverse of B / B in hybrid/mixed system
%         - type     : system type (hybrid or mixed)
%         - ip       : inner product name
%
% COMMENTS:
%  In the hybrid discretization, the matrices B, C and D appear as
%       [  B   C  D  ]
%       [  C'  O  O  ]
%       [  D'  O  O  ]
%
% SEE ALSO:
%   assembleWellSystem, solveIncompFlow, darcy, permTensor.

%{
#COPYRIGHT#
%}

% $Id: computeMimeticIP.m 3052 2009-10-22 10:42:37Z jrn $

error(nargchk(2, inf, nargin, 'struct'));

dim = size(G.nodes.coords, 2);
%--------------------------------------------------------------------------
%- Model and discretisation parameters ------------------------------------
%
perm = permTensor(rock, dim);
load Kreal;
Kreal     = reshape( Kreal, [1,G.cells.num]);
rock.perm = Kreal'*100*milli*darcy;
perm1     = permTensor(rock, dim);

[systemType, verbose, ip, ipname, ...
   computeIP, computeInverseIP] = setup(varargin{:});

computeDK = true ; %compute derivative of permeability
%--------------------------------------------------------------------------
%- Compute discrete gradient (C) and flux continuity (D) matrices ---------
%
cf = double(G.cellFaces(:,1));

if verbose,
   fprintf('Computing component matrices C and D ...\t')
   tic
end

cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
sgn    = 2*double(G.faces.neighbors(cf, 1) == cellno) - 1;

tocif(verbose)

%--------------------------------------------------------------------------
%- Derive geometric primitives per half-contact in given model ------------
%
areas      = G.faces.areas    (cf  );
centVecs   = G.faces.centroids(cf,:) - G.cells.centroids(cellno,:);
outNormals = bsxfun(@times, G.faces.normals(cf,:), sgn);
clear cf cellno sgn

%--------------------------------------------------------------------------
%- Preallocate SPARSE storage for BI, one value per half-contact ----------
%
dimProd = double(diff(G.cells.facePos));
sizVec  = sum(dimProd .^ 2);
ind1    = zeros([sizVec, 1]);
ind2    = zeros([sizVec, 1]);

if computeIP,        val_ip    = zeros([sizVec, 1]); end
if computeInverseIP, val_invip = zeros([sizVec, 1]); end
if computeDK,        val_dk    = zeros([sizVec, 1]); end %derivative permeability

%--------------------------------------------------------------------------
%- Compute half-contact numbering for each cell (-> SPARSE assembly) ------
%
[I1, I2] = localFaceIndices(dimProd);

% Current index to row/column in BI (or B), and current index to value.
ix = 0; ix2 = 0;

%--------------------------------------------------------------------------
% Main loop -- Compute BI and/or B for each cell in model -----------------
%

if verbose,
   fprintf('Computing cell inner products ...\t\t')
   tic
end
for k = 1 : G.cells.num,
   nF  = dimProd(k);
   pR  = ix  + 1 : ix  + nF;
   pR2 = ix2 + 1 : ix2 + nF^2;

   a = areas(pR);
   v = G.cells.volumes(k);

   K = reshape(perm(k,:), [dim, dim]);
   K1 = reshape(perm1(k,:), [dim, dim]);
   C = centVecs  (pR,:);
   N = outNormals(pR,:);

   % Compute inner product if caller specified 'mixed' method.
   if computeIP,
      val_ip(pR2)    = ip(a, v, K, C, N, false, false);
   end

   % Compute inverse inner product if caller specified 'hybrid' method.
   if computeInverseIP,
      val_invip(pR2) = ip(a, v, K, C, N, true, false);
   end
   
   % compute derivative of permeability
   if computeDK,
       val_dk(pR2) = ip(a, v, K1, C, N, false, true);
   end

   ind1(pR2) = I1{nF} + ix;
   ind2(pR2) = I2{nF} + ix;

   ix  = ix  + nF;
   ix2 = ix2 + nF.^2;
end
tocif(verbose)

%--------------------------------------------------------------------------
%- Assemble final block diagonal BI and/or B for complete model -----------
%
clear areas dimProd centVecs outNormals perm
if verbose,
   fprintf('Assembling global inner product matrix ...\t')
   tic
end
if computeIP,
   n    = size(G.cellFaces, 1);
   S.B  = sparse(ind1, ind2, val_ip,    n, n);
end
if computeInverseIP,
   n    = size(G.cellFaces, 1);
   S.BI = sparse(ind1, ind2, val_invip, n, n);
end
if computeDK, % compute derivative of permeability
   S.DK = sparse(ind1, ind2, val_dk,    n, n);
end
tocif(verbose)

if verbose && computeIP && computeInverseIP,
   fprintf('Max error in inverse = %g\n', ...
           norm(S.BI*S.B - speye(size(G.cellFaces,1)), inf));
end

%--------------------------------------------------------------------------
%- Define meta data on pressure system component matrices -----------------
%
S.ip    = ipname;
S.type  = systemType;

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------


%% Inverse of innerproduct - functions

function W = ip_simple(a, v, K, C, N, computeInverseIP, computeDK)
% These are not mutually inverse for general cells
t = 6/size(K,2) .* sum(diag(K));
if computeInverseIP,
   Q  = orth(bsxfun(@times, C, a));
   U  = eye(length(a)) - Q*Q';
   d  = diag(a);
   W  = (N*K*N' + t*(d * U * d)) ./ v;
else
   Q  = orth(bsxfun(@rdivide, N, a));
   U  = eye(length(a)) - Q*Q';
   di = diag(1 ./ a);
   W  = (C * (K \ C'))./v + (v / t)*(di * U * di);
end

%--------------------------------------------------------------------------

function W = ip_tpf(a, v, K, C, N, computeInverseIP, computeDK) %#ok
td = sum(C .* (N*K), 2) ./ sum(C.*C, 2);
if computeInverseIP,
   W = diag(td);
else
   W = diag(1 ./ td);
end
if computeDK  % derivative of permeability
    clear W;
    dk = sum(C.*C, 2) ./ sum(C .* (N*K), 2);
    W  = diag(dk);
    

end

%--------------------------------------------------------------------------

function W = ip_rt(a, v, K, C, N, computeInverseIP, computeDK) %#ok
Ki = inv(K);
Q  = kron(eye(size(K,2)), [1, 1]) ./ sqrt(2);
u  = diag(diag(Ki)) .* (v / 6);
di = diag(1 ./ a);
W  = C*(Ki*C')./v + di*Q'*u*Q*di;

if computeInverseIP, W = inv(W); end

%--------------------------------------------------------------------------

function W = ip_quasitpf(a, v, K, C, N, computeInverseIP, computeDK)
Q  = orth(C);
U  = eye(length(a)) - Q*Q';
W1 = N * K * N';

% ?? 2D ??
W  = (W1 + 2*U*diag(diag(W1))*U) ./ v;

if ~computeInverseIP, W = inv(W); end

%--------------------------------------------------------------------------

%% Setup helpers.

function s = id(s)
s = ['computeMimeticIP:', s];

%--------------------------------------------------------------------------

function [I1, I2] = localFaceIndices(nF)
found = reshape(find(accumarray(nF, 1)), 1, []);  % Row for FOR-loop below

I1 = cell([found(end), 1]);
I2 = cell([found(end), 1]);

for nf = found,
   I1{nf} =         repmat((1 : nf) .', [nf, 1]);
   I2{nf} = reshape(repmat( 1 : nf    , [nf, 1]), [], 1);
end

%--------------------------------------------------------------------------

function [Type, verbose, ipfun, ipName, IP, invIP] = setup(varargin)

opt = struct('InnerProduct', 'ip_simple', ...
             'Verbose',      false,       ...
             'Type',         'hybrid');
opt = merge_options(opt, varargin{:});

ipName  = opt.InnerProduct;
verbose = opt.Verbose;
Type    = opt.Type;

switch lower(Type),
   case 'hybrid',         [invIP, IP] = deal(true , false);
   case {'mixed','tpfa'}, [invIP, IP] = deal(false, true );
   case 'comp_hybrid',    [invIP, IP] = deal(true , true );
   otherwise,
      error(id('SystemType:Unknown'), ...
            'Unkown system type ''%s''.', Type);
end

ipNames = {'ip_simple', 'ip_tpf', 'ip_quasitpf', 'ip_rt'};
if isempty(find(strcmp(ipName, ipNames), 1)),
   error(id('InnerProduct:Unknown'), ...
         'Unknown inner product ''%s''.', ipName);
end
if verbose && strcmp(ipName, 'ip_rt'),
   disp(['NOTE: The mixed MFEM inner product is only valid for', ...
         ' Cartesian grids as given by cartGrid.'])
end

dispif(verbose, 'Using inner product: ''%s''.\n', ipName);
ipfun = str2func(ipName);
