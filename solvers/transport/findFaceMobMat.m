function [A_w, A_o] = ...
      findFaceMobMat(u_inx, mob, g_vec, v_darcy, neighbors, ...
                     normals, fluid)
%Initialize upwind saturation index matrices for face mobility.
%Uses known cell mobilities and output from function initFaceMob.
%
% SYNOPSIS:
%   [A_w, A_o] = findFaceMobMat(u_inx, mob, g_vec, v_darcy, ...
%                               neighbors, normals, fluid)
%
%
% PARAMETERS:
%   u_inx     - Index to cells for initial upwind face mobility not
%               considering gravity effects. Size = (#internal faces).
%
%   mob       - Cell mobility. (mob = fluid.mob(resSol)).
%
%   v_darcy   - Darcy flux for internal faces.
%
%   neighbors - Neighbors for internal faces.
%
%   normals   - Vector of face normals for internal faces.
%
%   g_vec     - Gravity contributions = abs[K(rho_1-rho_2)*g*n_z], for
%               each face, where n_z is the z-component of the normal.
%
%   fluid     - Fluid object as defined by function 'initFluid'.
%
% RETURNS:
%   A_w       - Index matrix for converting from cell mobility to face
%               mobility for water-phase (or primary phase).
%               facemob_w = A_w * mob.
%
%   A_o       - Index matrix for converting from cell mobility to face
%               mobility for oil-phase (or secondary phase).
%
% SEE ALSO:
%   initFaceMob, twophaseUpwFEGrav, twophaseUpwBEGrav.

%{
#COPYRIGHT#
%}

% $Id: findFaceMobMat.m 2962 2009-10-08 11:08:48Z ilig $


% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%   Given to cells i, j sharing a face e_ij we need to decide whether to
%   use the saturation from cell i (s_i) or cell j (s_j) when calculating
%   the mobility for e_ij. Let v_ij be the darcy flux over e_ij from cell i
%   to cell j. For cases without gravity, upstream weighting means that we
%   define the saturation of e_ij as
%
%           s(e_ij) = { s_i   if v_ij >= 0,
%                     { s_j   if v_ij <  0.
%
%   However, for a two-phase problem with gravity, the phase fluxes v_1 and
%   v_2 need not have the same direction as the Darcy flux: either v_1 or
%   v_2 can have oposite direction. Thus, in order to determine the
%   upstream mobility for each edge we need to know sign(v_1) and
%   sign(v_2).
%
%   Let f_k be the fractional flow of phase k and let
%      f_k(e_ij) := f_k(s_k(e_ij)), and  mob_k(s_k(e_ij)) := mob_k(e_ij).
%   Furhtermore let z be the depth and direction of gravity. The phase
%   fluxes for phases 1 and 2 over e_ij is given as
%      v_1(e_ij) = f_1(e_ij)*[v + K*(rho_1-rho_2)* mob_2(e_ij) *g*grad(z)]
%      v_2(e_ij) = f_2(e_ij)*[v - K*(rho_1-rho_2)* mob_1(e_ij) *g*grad(z)]
%
%   Observe that
%      sign(v_1) = sign[v + K*(rho_1-rho_2)* mob_2(e_ij) *g*grad(z)].
%
%   Since this equation involves mob_2(e_ij), it means that we need to know
%   the upstream face mobility for one of the phases in order to calculate
%   the sign of phase velocity for the other. We will explain how this is
%   possible below. First, let v_h be the phase velocity of the heaviest
%   phase and v_l the lightest phase (thus rho_h > rho_l), and let n_z be
%   the z-component of the normal pointing from cell i to j. We initialize
%   both phases with upstream weighting given by the Darcy velocity
%   (u_inx), and therefore need to change weighting for the faces where
%   sign(phase velocity)~=sign(Darcy velocity).
%
%   The possible scenarios for an edge e_ij are:
%
%   1. v >= 0 & n_z > 0  =>  v_h >= 0
%      sign(v_l) = sign(v - g_vec(e_ij)*mob_h(i)), if v_l < 0: change.
%   2. v <  0 & n_z > 0  =>  v_l <  0
%      sign(v_h) = sign(v + g_vec(e_ij)*mob_l(j)), if v_h > 0: change.
%   3. v >= 0 & n_z < 0  =>  v_l >= 0
%      sign(v_h) = sign(v + g_vec(e_ij)*mob_l(i)), if v_h < 0: change.
%   4. v <  0 & n_z < 0  =>  v_h <  0
%      sign(v_l) = sign(v - g_vec(e_ij)*mob_h(j)), if v_l > 0: change.
%
%   Notice that faces where n_z == 0 are not included since g_vec == 0 for
%   such faces. The result of the alogrithm is two arrays of cell-numbers,
%   h_inx and l_inx, giving the index to the upwind weighting for each
%   face. These arrays are used to build the upwind saturation index
%   matrices A_w and A_o.

   faceMob = mob(u_inx, :);

   h_inx = u_inx;
   l_inx = u_inx;
   
   % Recognize the heaviest phase:
   if fluid.rho(1)-fluid.rho(2) > 0 % 1. phase is heaviest
      h = 1; l = 2;
   else % 2. phase is heaviest
      h = 2; l = 1;
   end

   if nnz(g_vec),

      pos_v = ~(v_darcy < 0);  % v_darcy >= 0
      pos_n =   normals > 0;
      neg_n = ~pos_n;

      if norm(gravity()) > 0,
         inx1 =  pos_v & pos_n;
         inx2 = ~pos_v & pos_n;
         inx3 =  pos_v & neg_n;
         inx4 = ~pos_v & neg_n;
      else
         inx1 =  pos_v & neg_n;
         inx2 = ~pos_v & neg_n;
         inx3 =  pos_v & pos_n;
         inx4 = ~pos_v & pos_n;
      end

      
      %Four possible scenarios for a face:
      %1: v >= 0, n > 0
      if any(inx1)
         faces = find(inx1);
         cells2 = neighbors(inx1,2);
         a = (v_darcy(inx1) - faceMob(inx1,h).*g_vec(inx1)) < 0;
         l_inx(faces(a)) = cells2(a);
      end

      %2: v < 0, n > 0
      if any(inx2)
         faces = find(inx2);
         cells1 = neighbors(inx2,1);
         a = (v_darcy(inx2) + faceMob(inx2,l).*g_vec(inx2)) > 0;
         h_inx(faces(a))= cells1(a);
      end

      %3: v >= 0, n < 0
      if any(inx3)
         faces = find(inx3);
         cells2 = neighbors(inx3,2);
         a = (v_darcy(inx3) + faceMob(inx3,l).*g_vec(inx3)) < 0;
         h_inx(faces(a)) = cells2(a);
      end

      %4: v < 0 , n < 0
      if any(inx4)
         faces = find(inx4);
         cells1 = neighbors(inx4,1);
         a = (v_darcy(inx4) - faceMob(inx4,h).*g_vec(inx4)) > 0;
         l_inx(faces(a)) = cells1(a);
      end

   end %end if g_vec

   if h == 1 % Water (or the primary phase) is the heaviest phase.
      w_inx = h_inx;
      o_inx = l_inx;
   else
      w_inx = l_inx;
      o_inx = h_inx;
   end

   % Initialize upwind saturation index matrices
   n = numel(w_inx);
   A_w = sparse(1:n, double(w_inx), 1, n, size(mob,1));

   if o_inx == w_inx
      A_o = A_w;
   else
      A_o = sparse(1:n, double(o_inx), 1, n, size(mob,1));
   end
end
