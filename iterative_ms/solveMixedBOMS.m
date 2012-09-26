function [resSol, wellSol, incNorms, resNormsF, it] = ...
    solveMixedBOMS(resSol, wellSol, G, CG, part, rock, S, CS, fluid, ...
                   p0, dt, LS, I, J, Phi, Psi, varargin)

%--------------------------------------------------------------------------
%

opt = struct('LinSolve', @mldivide, 'tol', 1e-4, ... 
             'Verbose', true, 'maxIt', 100);
opt = merge_options(opt, varargin{:});

if opt.Verbose
    fprintf('   Solving MS Pressure Equation: \n');
end

[resFp,resF,resC] = deal(inf);
rs0     = resSol;
tolBOMS = opt.tol;


it = 0;
while resC > opt.tol && ...
      it < opt.maxIt %&& (isinf(resF) || (resF<0.99*resFp))
    [resSol, wellSol, incNorms, resNormsF, resNormsC, its] = ...
        iterateMixedBOMS(resSol, wellSol, G, ...
             CG, part, rock, S, CS, fluid, p0, dt, LS, ...
             I, J, Phi, Psi, 'tol', tolBOMS);
  
    inc = max(incNorms.flux, incNorms.pressure);
    resFp = resF;
    resF = max(resNormsF.f, resNormsF.g); 
    resC = max(resNormsC.f, resNormsC.g);
    it = it + its;
    if opt.Verbose
        fprintf('   %0*d: Relative Norms: Inc = %10.4e  ResF = %10.4e ResC = %10.4e\n', ...
                3, it, inc, resF, resC);
    end
    if inc < opt.tol, tolBOMS = 0.2*tolBOMS; end
end

incNorms.flux     = norm(resSol.faceFlux-rs0.faceFlux)/norm(resSol.faceFlux);
incNorms.pressure = norm(resSol.cellPressure-rs0.cellPressure)/(max(resSol.cellPressure)-min(resSol.cellPressure));