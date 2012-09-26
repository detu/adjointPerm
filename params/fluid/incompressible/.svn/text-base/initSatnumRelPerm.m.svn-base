function relperm = initSatnumRelPerm(T, varargin)
% Construct two-phase relperm evaluation function.
%
% SYNOPSIS:
%   relperm = initSatnumRelPerm(table)
%
% PARAMETERS:
%   table   - PVT/Relperm table structure as defined by function
%             'readpvt'.
%
% RETURNS:
%   relperm - Function for evaluating relative permeability curves
%             specified by SWOF (and SGOF).  Specifically, the call
%
%                [kr, dkr] = relperm(s)
%
%             evaluates the relative permeabilities (kr) and the
%             differentiated relative permeabilities (dkr) of the current
%             fluid composition (s, saturation).  
%
% REMARK:
%   adapted from initBlackOilRelPerm
%
% SEE ALSO:
%   initSatnumFluid, initBlackOilRelPerm, initBlackoilFluid, swof, sgof.
% 
%

%{
#COPYRIGHT#
%}

   opt = struct('verbose', false);	
   opt = merge_options(opt, varargin{:});
      
  	% Make struct of relperm tables
   if isfield(T, 'swof'),
      numSatnum = size(T.swof,2);
      kr_vec = cell(numSatnum, 2);    
      
      for i = 1:numSatnum   
         [krw, krow] = swof(struct('swof', T.swof{i}));
         kr_vec{i, 1} = krw;
         kr_vec{i, 2} = krow;     
      end      
   end
   
   % Return relperm according to satnum of the cells
   function [kr, varargout] = twophaseRelPerm(sol)      
      kr = zeros([size(sol.s,1), 2]);
      if nargout > 1,
         dkr = zeros([size(sol.s,1), 2]);
      end
      
      % default value e.g for plotting or computation of advection step
      if ~isfield(sol, 'satnum')
         if opt.verbose, disp('Using default value satnum = 1'), end;
         sol.satnum = true(size(sol.s));
         sol.satnumActive = 1;
      end
            
      for j = 1:numel(sol.satnumActive)
         k = sol.satnumActive(j);
                  
         [krw_s{1 : nargout}] = kr_vec{k, 1}(sol.s(sol.satnum(:,j))); 
         [kro_s{1 : nargout}] = kr_vec{k, 2}(1-sol.s(sol.satnum(:,j))); 
         kr(sol.satnum(:,j), :) = [krw_s{1}, kro_s{1}]; 
        
         if nargout > 1,
            dkr(sol.satnum(:,j), 1) = krw_s{2};
            dkr(sol.satnum(:,j), 2) = kro_s{2};  
         end
      end
      if nargout > 1,  varargout{1} = dkr; end
   end
  relperm = @twophaseRelPerm;
end

