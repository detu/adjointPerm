function [F, Jac] = twophaseJacobian(G, resSol, wellSol, rock, fluid, varargin)
%Implicit single point upwind solver for two-phase flow, including gravity.
%
% SYNOPSIS:
%   resSol = twophaseJacobian(G, resSol, rock, fluid)
%   resSol = twophaseJacobian(G, resSol, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwBEGrav solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time.  The nonlinear system of equations that must be
%   solved to move the solution from time=0 to time=tf, are solved using a
%   Newton-Raphson algorithm with line search to increase robustness.
%
%   In the case of failing to compute a solution using only a single step
%   over the interval [0,tf], an alternative strategy involving sub-steps
%   and step size control is employed.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   rock    - 
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%
% RETURNS:
%   F       - 
%
%   Jac     - 
%
% SEE ALSO:
%   implicitTransport, explicitTransport.

% TODO:
%   - implement gravity effects for pressure boundary and wells

%{
#COPYRIGHT#
%}

% $Date: 2009-10-14 08:05:00 +0200 (on, 14 okt 2009) $
% $Revision: 2994 $

   opt = struct('verbose', false, 'wells', [], 'src', [], 'bc', [], ...
                'use_fixed_directions', false);
   opt = merge_options(opt, varargin{:});

   assert (size(resSol.s,2)<3 || all(resSol.s(:,3)==0));
   
   
   % f, g and q are bound to nested functions Residual and Jacobian
   q              = computeTransportSourceTerm(resSol, wellSol, G, opt.wells, opt.src, opt.bc);
   [dflux, gflux] = getFlux(G, rock, fluid, resSol.faceFlux);
   
   constData = computeConstData(G, dflux, gflux);
   findMobIx = @(varargin) findFaceMobIx(constData, varargin{:});
   
   pv        = poreVolume(G, rock);
   q         = q(:);
   
   % System F(s) = 0 of non-linear equations (and the system's Jacobian
   % matrix) defining saturation equilibrium at a given time step.
   %
   % Here, we use
   %
   %               mw           krw(s)         kro(s)
   %    fw(s) = --------,  mw = ------,  mo = -------- 
   %             mw + mo          µw            µo
   %
   % With the matrices flux, grav, A_w and A_o, we can write an upwind
   % backward Euler discretisation of the Buckley-Leverett model as
   %
   %     s^(n+1) = s^n - (dt./pv)*(H(s^n+1) - max(q,0) - min(q,0)*f(s^n+1))
   %
   % where H(s) = (flux + grav*diag(A_o*lam_o(s))
   %                 *A_w*lam_w(s)./(A_w*lam_w(s)+A_o*lam_o(s)),
   %
   % lam_l is the mobility for phase l, f is Buckely-Leverett fractional
   % flow function, while A_o and A_w are index matrices that determine
   % the upstream mobility.
   %
   % The target function is
   %
   %    F(S) = S-S^{n-1} + (dt./pv)*(H(S) - max(q, 0)- min(q,0)*f(S))
   %
   % and the Jacobian is
   %
   %   dF(S) =  I + (dt./pv)(H'(S)  - min(q,0)*df(S))
      
   mob      = fluid.mob(resSol);
   [iw_fixed, io_fixed] = findMobIx(mob);
   
   
   
   %function setDirections(resSol)
   %   mob = fluid.mob(resSol);
   %   [iw_fixed, io_fixed] = findMobIx(mob);
   %end
   
   
   
   function F = Residual (resSol, s0, dt)
     %
     %  WARNING
     %  This function use dflux, gflux, q, pv, fluid and G from enclosing 
     %  function.
     %
     %  DISCRETISATION
     %
     %  Compute  F given by
     %
     %   F(s, s0, dt) = s-s0 + dt/pv\sum[H(s) - max(q,0) - min(q,0)*fw,
     %
     %   H_k(s) = \sum_i [fw_i(dflux_i + mo_i*gflux_i)]  (faces i of cell k)
     %
     %  where fw = mw/(mw+mo), mo and mw are mobilities, dflux and gflux
     %  are Darcy flux and gravity flux g(rho1-rho2)*face_normal,
     %  respectively.
     %
     %  Advancing an implicit upwind mobility weighted scheme one time step
     %  amounts to solving F(s, s0, dt) = 0, given s0 and dt.  
     
     %  The upwind mobility weighted flux is simply computed by evaluating 
     %  the water mobility mw in the upwind cell wrt water flux 
     %
     %    wflux = fw(dflux + mo*gflux).  
     %
     %  Similarily, we evaluate mo in the upwind cell wrt oil flux
     %
     %    oflux = fo(dflux - mw*gflux).  
     %
     %  The upwind directions for each phase are found using findMobIx, 
     %  and computeConstData.  The gravity flux is computed in getFlux.
     
      m = fluid.mob(resSol);
      if ~opt.use_fixed_directions,
         [iw, io] = findMobIx(m);
      else
         [iw, io] = deal(iw_fixed, io_fixed);
      end
      
      if ~any(gflux), assert(all(iw == io));end
      
      i   = all(G.faces.neighbors~=0, 2);
      ic  = reshape(double(G.faces.neighbors(i,:)), [], 1);
      sgn = [ones(sum(i), 1); -ones(sum(i),1)];
      io  = reshape(repmat(double(io), [1, 2]), [], 1);
      iw  = reshape(repmat(double(iw), [1, 2]), [], 1);
      
      
      fw_cell = m(:,1) ./ sum(m,2);
      
      m_face = [m(iw, 1) m(io,2)];
      fw_face = m_face(:,1)./sum(m_face,2);
      
      F = resSol.s(:,1) - s0;
      F = F + dt.*(1./pv).*(...
          -(max(q,0) + min(q,0).*fw_cell)...
          + accumarray( ic, sgn.*fw_face.*(dflux+gflux.*(m_face(:,2))) )...
                           );
      
   end

   function J = Jacobian (resSol, dt)
   %
   %  This function use dflux, gflux, q, pv, fluid and G from
   %  enclosing function.
   %
   %
   %  At an interface, we compute the Jacobain of F as follows,
   %
   %    dF      dt   mo     1
   %    --- =   --  -----  -----  (dflux + mo*gflux)
   %    dmw     pv  mw+mo  mw+mo
   %
   %            dt  fo
   %        =   --  --  (dflux + mo*gflux)
   %            pv  mt
   %
   %
   %    dF      dt  fw
   %    --- = - --  --  (dflux - mw*gflux)
   %    dmo     pv  mt
   %
   %
   %    dF    dF  dmw    dF  dmo
   %    --  = --- ---  + --- ---
   %    dS    dmw dS     dmo dS
   %
   %  If the mobilities mw and mo are evaluated in grid cells iw and io,
   %  respectively, then the respective (matrix) contributions of the each
   %  term is added to elements (i, iw) and (i, io) of the Jacobian matrix.
   %
   %
   %
      
      
      
      mob = fluid.mob(resSol);
      if ~opt.use_fixed_directions,         
         [iw, io] = findMobIx(mob);
      else
         [iw, io] = deal(iw_fixed, io_fixed);
      end
      
      i   = all(G.faces.neighbors~=0, 2);
      id = (1:G.cells.num)';      
      
      if ~any(gflux), 
         assert(all(iw == io));
         %assert(all(iw == G.faces.neighbors(i,1)));
      end
     
      sgn = [ones(sum(i), 1); -ones(sum(i),1)];
      
      ic  = reshape(double(G.faces.neighbors(i,:)), [], 1);      
      io  = reshape(repmat(double(io), [1, 2]), [], 1);
      iw  = reshape(repmat(double(iw), [1, 2]), [], 1);
      
            
      % [(d/ds) lam_w (d/ds) lam_o]
      dmob = bsxfun(@rdivide, fluid.dkr(resSol),fluid.mu);

      % On face
      mobf = [mob(iw, 1) mob(io,2)];
      ff   = bsxfun(@rdivide, mobf, sum(mobf,2));
      dw_fo    = ff(:,2).*dmob(iw,1);
      do_fw    = ff(:,1).*dmob(io,2);

      % In cell
      f   = bsxfun(@rdivide, mob, sum(mob,2));
      df  = (f(:,2).*dmob(:,1) - f(:,1).*dmob(:,2))./sum(mob,2);
      
      
      d = 1-dt./pv.*(min(q, 0).*df);
      dFdSiw = dt./pv(ic).*( sgn.*(dflux + gflux.*mobf(:,2)).*dw_fo./sum(mobf,2));
      %dFdSio = dt./pv(ic).*(-sgn.*(dflux + gflux.*mob(io,2)).*do_fw./sum(mobf,2) + sgn.*gflux.*do_fw);
      dFdSio = dt./pv(ic).*(-sgn.*(dflux - gflux.*mobf(:,1)).*do_fw./sum(mobf,2));
      
    
      
      J = sparse([id; ic; ic], [id; iw; io], [d; dFdSiw; dFdSio], ...
                 G.cells.num, G.cells.num);
      
    
      
      %{
      a = iw == ic; 
      b = io == ic;
      d = d + accumarray([ic(a);     ic(b)], ...
                         [dFdSiw(a); dFdSio(b)], [G.cells.num, 1]);
              
      J = sparse([id; ic(~a); ic(~b)], [id; iw(~a); io(~b)], [d; dFdSiw(~a); dFdSio(~b)], ...
                 G.cells.num, G.cells.num);
      %}
   end
   
   % Return function handles
   F   = @Residual;
   Jac = @Jacobian;
end

%--------------------------------------------------------------------------


function [f, g] = getFlux(G, rock, fluid, darcyFaceFlux)
% return total velocity and gravity face flux.
   is_int = all(G.faces.neighbors~=0, 2);

   [nc, nf] = deal(G.cells.num, sum(double(is_int)));
   cellNo   = rldecode(1 : nc, double(G.cells.numFaces), 2) .';

   % Indices to (internal) half-faces.
   cIntFInx = is_int(G.cellFaces(:,1));

   % Indices of internal face corresponding to cIntFInx.
   globfix  = G.cellFaces(cIntFInx, 1);

   % rho(1) - rho(2) is always correct (even when rho(2) > rho(1)), because
   % during transport (using, e.g., the 'twophaseUpwBEGrav' transport
   % solver) we only modify the *first* saturation component.
   g   = gravity() * (fluid.rho(1) - fluid.rho(2));
   dim = size(G.nodes.coords, 2);

   harm          = zeros([G.faces.num, 1]);
   renum         = zeros([G.faces.num, 1]);
   renum(is_int) = 1 : nf;

   % nKg == n' * K * g for all cellfaces.
   [K, r, c] = permTensor(rock, dim);
   nKg = sum(G.faces.normals(G.cellFaces(:,1), r) .* ...
             bsxfun(@times, K(cellNo,:), g(c)), 2);

   % Compute harmonic average of nKg on all *internal* faces.
   harm(is_int) = 1 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));

   f = darcyFaceFlux(is_int);
   f = [f;f];
   
   g = harm(is_int);
   g = [g;g];
end

%--------------------------------------------------------------------------

function constData = computeConstData(G, myff, mygg)  
%
% Given fixed dflux and gflux, we can compute upwind directions for 
% phase w given that dflux and gflux have the same sign.  If dflux and 
% gflux have opposite signs, then the upwind directions for phase o can be
% computed.
%
% We store these precomputed values to speed up calculations of upwind
% direcions in findFaceMobIx.
   intern     = all(G.faces.neighbors ~= 0, 2);
   v_darcy    = myff(1:end/2);
   g_vec      = mygg(1:end/2);
   ineighbors = G.faces.neighbors(intern,:);
   constData.ineighbors = ineighbors;
   
   
   vec = v_darcy .* g_vec > 0;
    
   iw = zeros(sum(intern), 1);
   io = iw;
   
   ineighbors(v_darcy<0, :) = ineighbors(v_darcy<0, [2,1]);
   iw( vec) = ineighbors( vec, 1);
   io(~vec) = ineighbors(~vec, 1);
   
   
   constData.v_darcy = v_darcy;
   constData.g_vec   = g_vec;
   constData.vec     = vec;
   constData.io      = io;
   constData.iw      = iw;
   
   
end

                       
%--------------------------------------------------------------------------

function [IW, IO] = findFaceMobIx(constData, mob)
%Compute upwind cell index for each phase flux (internal to grid). 
%
% SYNOPSIS:
%   [iw, io] = findFaceMobIx(constData, mob)
%
%
% PARAMETERS:
%   constData - struct computed in computeConstData
%
%   mob       - Cell mobility. (mob = fluid.mob(resSol)).
%
% RETURNS:
%   iw       - Index for converting from cell mobility to face
%               mobility for water-phase (or primary phase).
%               facemob_w = mob(iw,1).
%
%   io       - Index for converting from cell mobility to face
%               mobility for oil-phase (or secondary phase).
%
% SEE ALSO:
%   initFaceMob, twophaseJacobian.

%{
#COPYRIGHT#
%}

% $Id: twophaseJacobian.m 2994 2009-10-14 06:05:00Z jrn $


% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  For each pair of Darcy flux dflux and gravity flux gflux=g·n(rho1-rho2),
%  we can compute the phase velocities for water (w) and oil (o),
%
%           mw
%     vw = -----  (dflux + mo*gflux)
%          mw+mo
%
%
%           mo
%     vo = -----  (dflux - mw*gflux)
%          mw+mo
%
%   
% To use the upwind mobility flux approximation, we evaluate the water
% mobility mw in the upwind direction wrt vw, and the oil mobility mo in
% the upwind direction wrt vo.
%
% This is solved by noting that if dflux and gflux are positive(negative),
% the vw is positive(negative), regardless of the value of mo.  Thus, the
% upwind direction of water is known, mw may be evaluated, vo computed and
% the upwind direction for oil can be found.
%
% If dflux is positive(negative) and gflux is negative(positive), vo is
% positive(negative) and the rest may be computed.

   
   iw = constData.iw;
   io = constData.io;
   
   mw = zeros(size(iw));
   mo = mw;
   
   vec = constData.vec;
   mw( vec) = mob(iw( vec), 1);
   mo(~vec) = mob(io(~vec), 2);
   
   vw = constData.v_darcy +  constData.g_vec.*mo;
   vo = constData.v_darcy -  constData.g_vec.*mw;
   
   ineighbors = constData.ineighbors;
   ineighbors(vw<0, :) = ineighbors(vw<0, [2,1]);
   iw(~vec) = ineighbors(~vec, 1);
   
   ineighbors = constData.ineighbors;
   ineighbors(vo<0, :) = ineighbors(vo<0, [2,1]);
   io(vec) = ineighbors(vec, 1);


   IW = iw;
   IO = io;
end
