function [resSol, wellSol, incNorms, resNorms] = ...
   solveMixedBO(resSol, wellSol, G, rock, S, fluid, dt, varargin)

%--------------------------------------------------------------------------
%

opt = struct('bc', [], 'src', [], 'wells', [], ...
             'LinSolve', @mldivide, 'tol', 1e-4, ... 
             'Verbose', true, 'maxIt', 100);
opt = merge_options(opt, varargin{:});

tol = opt.tol;
maxIt = opt.maxIt;

if opt.Verbose
    fprintf('Solving Pressure Equation: \n');
end

LS = initMixedBO(G, opt.bc, opt.src, opt.wells);

faceFlux = [resSol.faceFlux; vertcat( wellSol.flux )];
p    = resSol.cellPressure;
lam  = [resSol.facePressure; vertcat( wellSol.pressure )];
lamN = lam( [LS.fluxFacesR; LS.fluxFacesW] );

it = 0;
[inc,res] = deal(inf);
p0 = resSol.cellPressure;
while max(inc,res) > tol && it < maxIt

   LS = setupMixedBO(resSol, wellSol, G, rock, S, fluid, p0, dt, LS);

   incNorms     = [];
   resNorms.f = norm(LS.f)/(max(p)-min(p));
   resNorms.g = norm(LS.g)/norm(faceFlux);
   if max(resNorms.f, resNorms.g) > opt.tol
      [flux_inc, p_inc, lamN_inc] = ...
         solveMixedLinSys(LS.B, LS.C, LS.DN, LS.P, LS.f, LS.g, LS.hN, LS.Do);
      LS.solved = true;

      faceFlux = faceFlux + flux_inc;
      p    = p + p_inc;
      lamN = lamN + lamN_inc;

      incNorms.flux     = norm(flux_inc)/norm(faceFlux);
      incNorms.pressure = norm(p_inc)/(max(p)-min(p));

      [resSol, wellSol] = ...
         packSol(G, opt.wells, resSol, wellSol, faceFlux, p, lamN, opt.bc);
      
      it = it + 1;
   else
      incNorms.flux     = 0;
      incNorms.pressure = 0;
      LS.solved = false;
   end

   inc = max(incNorms.flux, incNorms.pressure);
   res = max(resNorms.f, resNorms.g);
   if opt.Verbose
      fprintf('   %0*d: Relative Norms: Inc = %10.4e  Res = %10.4e\n', 3, it, inc, res);
      %fprintf('   %0*d: Max RelResNorm = %10.4e\n', 3, it, res);
   end
end

%--------------------------------------------------------------------------

function [resSol, wellSol] = packSol(G, W, resSol, wellSol, flux, ...
                                     p, lamN, bc)
fluxFacesR = getBCType(G, W, bc);
facePressure = zeros(G.faces.num, 1);
if ~isempty(bc)
    pf = strcmpi('pressure', bc.type');
    presFaces = double(bc.face(pf));
    facePressure(presFaces) = bc.value(pf);
end
facePressure(fluxFacesR) = lamN(1:nnz(fluxFacesR));

resSol.facePressure = facePressure;
resSol.cellPressure = p;
resSol.faceFlux     = flux(1:G.faces.num);
resSol.cellFlux     = faceFlux2cellFlux(G, flux);

%wells
r_inx = G.faces.num;
p_inx = numel(fluxFacesR);
for k = 1:numel(W)
    if strcmpi(W(k).type, 'pressure')
       wellSol(k).pressure = W(k).val;
    else
        wellSol(k).pressure = lamN(p_inx+1);
        p_inx = p_inx+1;
    end
    nwc = numel(W(k).cells);
    wellSol(k).flux = - flux(r_inx + (1:nwc));
    r_inx = r_inx + nwc;
end