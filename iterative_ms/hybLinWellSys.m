function [B, f, g, h] = ...
    hybLinWellSys(resSol, wellSol, G, W, fluid)
if isempty(W)
    [B,f,h] = deal([]);
    g = sparse(G.cells.num, 1);
else
    nperf = cellfun(@numel, { W.cells });
    wc    = vertcat(W.cells);   np = numel(wc);
    rates = vertcat( wellSol.flux );
    
    [c, rho, mu, u] = fluid.pvt(resSol.cellPressure(wc), resSol.z(wc));
    s      = bsxfun(@rdivide, u, sum(u,2));
    mob    = fluid.relperm(s(wc)) ./ mu;
    Lti    = 1./sum(mob, 2);
    rhoTot = sum(rho.*mob ,2).*Lti;
    a1     =  sum( c .* mob, 2) .* Lti;
    
    % Diagonal transmissibility matrix (inner product).
    B   = sparse(1:np, 1:np, Lti .* 1./( vertcat(W.WI) ), np, np);
    
    % Which (and what) are the prescribed well bottom-hole pressures?
    dF   = strcmpi('bhp', { W.type } .');
    
    % Form linsys rhs contributions, fW -> pressure, hW -> rate.
    f = -vertcat(W.val);  f(~dF) = 0;  % Remove rates
    h = -vertcat(W.val);  h( dF) = 0;  % Remove pressures
    g = sparse(G.cells.num, 1);
    g(wc) = a1 .* rates .* (B * rates);
    
    % Expand well pressure rhs to each perforation, adjust for gravity
    % effects if applicable (i.e., when NORM(gravity())>0 and W.dZ~=0).
    dp  = computeDp(W, rho, rhoTot, wellSol);
    f   = rldecode(f, nperf) - dp;
end

function dp = computeDp(W, rho, rhoTot, wellSol)
% Find pressure drop distribution along well. Assume all phases flow freely
% at infinite speed such that segment volume fraction of phase a is equal to
% inflow fraction of phase a. Pressure drop in a segment is then rho*g*Dz,
% where rho is the mixture density. Further for now make simplifying assumption
% that densities in segment equal densities in well-cell.
nperf = cellfun(@numel, { W.cells });
dp    = zeros(prod(nperf), 1);
grav  = gravity();
if norm(grav)>0
    inx = 0;
    for k = 1 : numel(W)
        wc = W(k).cells;
        numSeg   = numel(wc);
        inxs = inx + (1:nperf(k))';
        inx  = inx + nperf(k);
        if numSeg > 1
            resFlux  = wellSol(k).flux;
            totFlux  = sum( resFlux );
            if any(sign(resFlux) == -sign(totFlux)  )
                warning('Both in- and out-flux in same well.')
            end
            wellFlux = [totFlux; totFlux - cumsum( resFlux(1:end) )];
            
            rhoRes  = rhoTot(inxs);
            rhoInj  = rho(inxs(1),:) * W(k).Comp_i';
            
            res2segFlux    = resFlux .* (resFlux < 0);
            seg2segPosFlux = wellFlux .* (wellFlux > 0); %heel to toe (i -> i+1)
            seg2segNegFlux = wellFlux .* (wellFlux < 0); %toe to heel (i -> i-1)
            
            rhoSeg = rhoRes;
            
            % iterate until fixed
            rhoSegPrev = rhoSeg*0;
            while rhoSegPrev ~= rhoSeg
                rhoSegPrev = rhoSeg;
                rhoSeg = res2segFlux .* rhoRes + ...
                    seg2segPosFlux(1:numSeg) .* [rhoInj; rhoSeg(1:numSeg-1)] +...
                    seg2segNegFlux(2:numSeg+1) .* [rhoSeg(2:numSeg); 0];
            end
            perfDZ = W(k).dZ(2:numSeg) - W(k).dZ(1:numSeg-1);
            perfDp = norm(gravity()) * ...
                [0; .5 * perfDZ * (rhoSeg(1:numSeg-1) + rhoSeg(2:numSeg))];
            dp(inxs) = cumsum( perfDp );
        end
    end
end


            




