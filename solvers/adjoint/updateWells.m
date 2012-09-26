function [W] = updateWells(W, scheduleStep)
% Update wells based on schedule time step
numWells   = numel(W);
multiscale = isfield(W(1), 'CS');

for k = 1 : numWells
    type      = scheduleStep.types{k};
    val       = scheduleStep.values(k);
    W(k).type = type;
    W(k).val  = val;  
    if strcmp(type, 'bhp'),
      W(k).S.RHS.f = - val(ones([W(k).S.sizeB(1), 1]));
      W(k).S.RHS.h = 0;                % never used
      if multiscale
          numCoarseCells = length( W(k).coarseCells );
          W(k).CS.RHS.f = -W(k).val(ones([numCoarseCells, 1]));
          W(k).CS.RHS.h = 0;      % never used
      end 
   elseif strcmp(type, 'rate'),
      W(k).S.RHS.f = zeros([W(k).S.sizeB(1), 1]);
      W(k).S.RHS.h = - val;
      if multiscale
          numCoarseCells = length( W(k).coarseCells );
          W(k).CS.RHS.f = zeros([numCoarseCells, 1]);
          W(k).CS.RHS.h = -W(k).val;
      end 
   end
end
