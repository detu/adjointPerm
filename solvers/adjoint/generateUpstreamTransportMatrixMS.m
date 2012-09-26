function [A, qPluss, signQ] = ...
  generateUpstreamTransportMatrixMS(G, S, W, SG, resSol, wellSol, varargin)
% generateUpstreamTransportMatrixMS for use in saturation solver
%
% SYNOPSIS: 
%   [A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, ...
%                                           resSol, wellSol, pn1, pv1, ...)
%
% DESCRIPTION:
%   Generates sparse matrix A(cellFlux) which is used in transport solver 
%   s1 = s0 + dt*Dv*(Af(s) + q+). 
%   Assumes no-flow boundary conditions on all cell-faces.
% 
%   signQ is useful for differentation of max(q, 0)/min(q, 0) wrt to q,
%   when q is zero
%
% REQUIRED PARAMETERS: 
%  
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   Transpose     - if true, the transpose of a is given (default false)
%
%   VectorOutput  - if true, output A is a struct with fields 'i', 'j' and
%                   'qMinus', such that A = sparse(i, j, -cellFlux) +
%                   diag(qMinus) (NOTE MINUSES) Default value is false
%
%   RelativeThreshold  - Considers values below
%                        max(abs(cellflux))*RealtiveThreshold 
%                        as zero (default 0)
%

% ---------------------------------------------
opt = struct('Transpose'         , false, ...
             'VectorOutput'      , false, ...
             'RelativeThreshold' , 0);
opt = merge_options(opt, varargin{:});

% ---------------------------------------------

cellFlux    = resSol.cellFlux;
inx         = (cellFlux ~= 0); 
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Cw          = vertcat(wellSys.C);

cfns  = cellFaceNeighbors(G, inx);
ii    = cfns(:,1);
q     = Cw'*wellRates;

%if opt.RelativeThreshold > 0
%    maxV              = max( abs(cellFlux) );
%    belowThresholdInx = ( abs(cellFlux) < opt.RelativeThreshold/maxV );
%    cellFlux(belowThresholdInx)  = 0;
%    belowThresholdInx = ( abs(q) < opt.RelativeThreshold/maxV );
%    q(belowThresholdInx)         = 0;
%end

% only allow interior faces
negFlux     = and( cellFlux(inx) < 0, cfns(:,2));  

%extInd      = (ii == 0);
%ii(extInd)  = cfns(extInd, 2);

part   = SG.cells.subCells * (1 : SG.cells.num)';

jj          = ii;
jj(negFlux) = cfns(negFlux, 2);

qMinus        = q;
qMinus(q > 0) = 0;

% -----------------------------------

if opt.VectorOutput
    A.i        = part( ii ); 
    A.j        = part( jj );
    A.qMinus   = SG.cells.subCells'*qMinus;
    A.cellFlux = cellFlux( inx ); 
else
    numC  = SG.cells.num;
    if opt.Transpose
        A = sparse(part( jj ), part( ii ), -cellFlux(inx), numC, numC) ...
            + spdiags(SG.cells.subCells'*qMinus, 0, numC, numC);
    else
        A = sparse(part( ii ), part( jj ), -cellFlux(inx), numC, numC) ...
            + spdiags(SG.cells.subCells'*qMinus, 0, numC, numC);
    end
end

if nargout > 1
    qPluss = SG.cells.subCells'*(q - qMinus);
end

if nargout > 2
    Dw    = blkdiag( wellSys.D );
    if isfield(W, 'sign')        
        signs = vertcat(W.sign);
    else
        signs    = ones( numel(W), 1 );
        totRates = Dw'*wellRates;
        signs( totRates < 0 ) = -1;
    end
    signQ       = SG.cells.subCells' * (Cw'*Dw*signs);
end
        
            
function cfns = cellFaceNeighbors(G, inx)
% finds neighbor pairs according to cellface - cells
% [cellface-cells celface-cell-neigbors

cfns     = double( G.faces.neighbors(G.cellFaces(inx,1), :) );
cellNo   = rldecode(1:G.cells.num, double(G.cells.numFaces), 2)';
sgn      = 2*(G.faces.neighbors(G.cellFaces(:,1), 1) == cellNo)-1;
flipRowInx  = ( sgn(inx) < 0 );
cfns(flipRowInx, :) = cfns(flipRowInx, [2 1]);
extInd      = find(cfns(:,1)==0);
cfns(extInd, :) = cfns(extInd, [2 1]);
return


