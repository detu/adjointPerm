function [xr, xw, incNorms, resNormsF, resNormsC, it] =  ...
    iterateMixedBOMS(xr, xw, G, CG, part, rock, S, CS, fluid, ...
                     p0, dt, LS, I, J, Phi, Psi, varargin)


opt = struct('tol', 0, 'LinSolve', @mldivide);
opt = merge_options(opt, varargin{:});

it = 0;

LS = setupMixedBO(xr, xw, G, rock, S, fluid, p0, dt, LS);

faceFlux = [xr.faceFlux; vertcat(xw.flux)];
p = xr.cellPressure;

%------ check fine-scale residual -----------------------------------------
resNormsF.f = norm(LS.f)/(max(p)-min(p));
resNormsF.g = norm(LS.g)/norm(faceFlux);
if max(resNormsF.f,resNormsF.g)<opt.tol
    incNorms.flux     = 0;
    incNorms.pressure = 0;
    resNormsC.f       = -inf;
    resNormsC.g       = -inf;
    return
end;

%------ restriction operators --------------------------------------------- 
DM = mobScaling(G, CG, CS, part, LS.mobt);
Phi = DM * Phi;

neumGP = [LS.fluxFacesR; LS.fluxFacesW];

[A, b, inx] = coarseSys(LS, I, J, Phi, Psi, neumGP);
resNormsC.f = norm(b(inx{1}))/(max(p)-min(p));
resNormsC.g = norm(b(inx{2}))/norm(faceFlux);
% if max(resNormsC.f,resNormsC.g)<opt.tol
%     incNorms.flux     = 0;
%     incNorms.pressure = 0;
%     return
% end;

x = opt.LinSolve(A,full(b)); LS.solved = true;
it = 1; % only count iterations where we solve a system..

sol = { x(inx{1}), -x(inx{2}), x(inx{3}) };

% set non-Neumann pressures to zero !!!!!!!!!!!!!
xr.facePressure(~neumGP) = 0;
[xr_inc, xw_inc] = pack_solution(xr, xw, G, CG, part, CS, sol{:}, ...
                                Phi, Psi, neumGP, LS.wells);
xr = resPlus(xr, xr_inc);
xw = wellPlus(xw, xw_inc);

pw = strcmpi('pressure',[LS.bc.type])';
xr.facePressure(LS.bc.face(pw)) = LS.bc.value(pw);

incNorms.flux     = norm(xr_inc.faceFlux)/norm(xr.faceFlux);
incNorms.pressure = norm(xr_inc.cellPressure)/(max(xr.cellPressure)-min(xr.cellPressure));
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------
function resSol = resPlus( resSol, rs)
fn = {'cellPressure', 'facePressure', 'cellFlux', 'faceFlux'};
for k = 1:numel(fn)
    resSol.(fn{k}) = resSol.(fn{k}) + rs.(fn{k}); 
end
end

function wellSol = wellPlus( wellSol, ws)
for wNr = 1 : numel(wellSol)
    fn = fieldnames(wellSol(k));
    for k = 1:numel(fn)
        wellSol(wNr).(fn{k}) = wellSol(wNr).(fn{k}) + ws(wNr).(fn{k});
    end
end
end

function [A, b, inx] = coarseSys(LS, I, J, Phi, Psi, neum_f)
    % Symmetric version based on reduction
    %  / Psi  0  0 \
    % |  Phi  I  0  |
    %  \  0   0  J /
    % thus no basis for boundary pressure yet
    
    JN_f = J(neum_f, :);
    neum_c = logical( sum( abs(JN_f) ) );
    JN   = JN_f(:, neum_c);
    
    B  = Psi' * LS.B * Psi;% + Phi' * P * Phi;
    C  = Psi' * (LS.C * I);% - Phi' * P * I;
    Mt = I' * LS.P * Phi;
    DN = Psi' * LS.DN * JN;
    P  = I' * LS.P * I;
    
    f   = Psi' * LS.f;% - Phi' * g;
    g   = I' * LS.g;
    hN  = JN' * LS.hN;
    
    nf = numel(f); ng = numel(g); nh = numel(hN);
    
    A = [ B          C              DN
          C'-Mt      P         sparse(ng, nh)
          DN'  sparse(nh, ng)  sparse(nh, nh)];
    b = [f; g; hN];
    inx = { 1:nf, nf + (1:ng), nf + ng + (1:nh) };
end



function DM = mobScaling(G, CG, CS, part, mob_f)
    %cumbersome!
    mobt0_af = cellfun(@(x) x{6}, { CS.basisP{CS.activeFaces} }, ...
                  'UniformOutput', false);
    cells_af = cellfun(@(x) x{4}', { CS.basisP{CS.activeFaces} }, ...
                  'UniformOutput', false);
    mobt0    = zeros( CG.cells.num, 1 );
    mobt0( vertcat( cells_af{:} ) ) = vertcat( mobt0_af{:} ); %coarse cell mob
   
    mobt  = accumarray(part, mob_f.* G.cells.volumes, [CG.cells.num, 1]) ./ ...
            accumarray(part,        G.cells.volumes, [CG.cells.num, 1]);
    ms    = mobt0./mobt;
    
    DM = sparse(1:G.cells.num, 1:G.cells.num, ms(part));
end


%--------------------------------------------------------------------------

function [xr, xw] = pack_solution(xr, xw, G, CG, p, CS, flux, pres, ...
                                  lam, Phi, Psi, neumf, w)
   % [nsub, sub] = get_subfaces(G, CG, CS);
   actB = CS.activeCellFaces;
   actF = CS.activeFaces;

   %-----------------------------------------------------------------------
   % Package reservoir (coarse) solution structure ------------------------
   %
   %switch lower(opt.Solver),
   %   case 'hybrid',
   %      nhf      = numel(actB);
   %      flux_sgn = rldecode([1, -1], [nhf, size(Phi,2) - nhf], 2) .';
   % 
   %      xr.blockFlux = flux(1 : numel(actB));
   %      fluxW        = flux(numel(actB) + 1 : end);
   %   case 'mixed',
         nf       = numel(actF);
         flux_sgn = rldecode([1, -1], [nf, size(Phi,2) - nf], 2) .';

         orient       = get_orientation(CG, CS);
         renum        = zeros([CG.faces.num, 1]);
         renum(actF)  = 1 : numel(actF);

         xr.blockFlux = orient .* flux(renum(CG.cellFaces(actB,2)));
         fluxW        = flux(numel(actF) + 1 : end);
   %end

   cellNo = rldecode(1 : G.cells.num, double(G.cells.numFaces), 2).';

   xr.blockPressure   = pres;
   xr.cellPressure(:) = pres(p) + Phi * (flux .* flux_sgn);
   xr.facePressure(:) = accumarray(G.cellFaces(:,1), ...
                                   xr.cellPressure(cellNo)) ./ ...
                        accumarray(G.cellFaces(:,1), 1);
   xr.facePressure(~neumf) = 0;

   %pw = strcmpi('pressure',[opt.bc.type])';
   %xr.facePressure(opt.bc.face(pw)) = opt.bc.value(pw);
                 
   %if strcmpi(opt.Solver, 'hybrid'),
   %   % The 'mixed' solver does not generate reasonable coarse interface
   %   % pressure on internal coarse faces.  Therefore, only insert coarse
   %   % interface pressures computed by the solver in the 'hybrid' case.
   %   %
   %   xr.facePressure(sub) = rldecode(lam(1 : numel(actF)), nsub);
   %end

   xr.cellFlux(:)     = Psi * flux;
   xr.faceFlux(:)     = cellFlux2faceFlux(G, xr.cellFlux);

   %-----------------------------------------------------------------------
   % Recover fine scale well fluxes and pressures. ------------------------
   %
   %w    = opt.wells;
   lamW = lam(numel(actF) + 1 : end);  assert (numel(lamW) == numel(w));

   nw = numel(w);
   i  = 0;
   for k = 1 : nw,
      rates = w(k).CS.rates;

      xw(k).flux(:)  = - full(rates * fluxW(i + 1 : i + size(rates,2)));
      xw(k).pressure = lamW(k);

      i = i + size(rates,2);
   end
end

%--------------------------------------------------------------------------

function orient = get_orientation(cg, cs)
   c1     = cg.faces.neighbors(cg.cellFaces(:,2), 1);
   orient = 2*double((cg.cellFaces(:,1) == c1)) - 1;
   orient = orient(cs.activeCellFaces);
end

%--------------------------------------------------------------------------
