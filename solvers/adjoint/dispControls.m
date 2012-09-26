function [] = dispControls(controls, schedule, varargin)


fprintf('\n----------------- DISPLAYING CONTROL VARIABLES ----------------\n')

numC    = numel(controls.well);
cw      = [controls.well(:).wellNum];
names   = [schedule(1).names(cw)];
types   = {controls.well(:).type};
minMax  = vertcat(controls.well.minMax);

fprintf('%9s%9s%9s%15s\n', 'Var', 'Name', 'Type', 'MaxMin')
for k = 1: numC
    var     = sprintf('u_%d', k);
    vars{k} = var;
    fprintf('%9s%9s%9s%15s\n', var, names{k}, types{k}, mat2str(minMax(k, :)) );
end

ec = controls.linEqConst;
if ~isempty(ec)
    fprintf('\nLinear equality constraints: \n')
    numEC = size(ec.A, 1);
    for k = 1: numEC
        A = ec.A(k, :); b = ec.b(k);
        [r, c] = size(A);
        lstart = false;
        for k1 = 1:c
            if A(k1) ~= 0
                if lstart, fprintf(' + '); end
                if A(k1) ~= 1, fprintf(num2str(A(k1))); end
                fprintf(vars{k1});
                lstart = true;
            end
        end
        fprintf(' = ');
        fprintf(num2str(b));
        fprintf('\n');
    end
end




