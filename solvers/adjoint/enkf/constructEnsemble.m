%----------------------------------------------------------------------------------
% SYNOPSIS:
%
% DESCRIPTION:
%   Construct true state vector xt and ensemble of states A
%
% PARAMETERS:
%
% RETURNS:
%   xt, A and nstat
%
% AUTHOR: 
%   Olwijn Leeuwenburgh (TNO)
%
% Copyright 2012 TNO. This file is part of the EnKF module. EnKF is free code.
%---------------------------------------------------------------------------------- 
ns = []; is = 0;

% permeability
if ~isempty(strmatch('prm',variables))
    prm = rock.perm(:,1);
    if ne > 0
        prme = zeros(numel(prm),ne);
        for j=1:ne
            prme(:,j) = rockE{j}.perm(:,1);
        end
    end
    is = is + 1;
    ns(is) = numel(prm);    
end

% porosity
if ~isempty(strmatch('por',variables))
    por = rock.poro;
    if ne > 0
        pore = zeros(numel(por),ne);
        for j=1:ne
            pore(:,j) = rockE{j}.poro; 
        end
    end
    is = is + 1;
    ns(is) = numel(por);
end

% saturation
if ~isempty(strmatch('sat',variables))
    sat = rSol.s; 
    if ne > 0
        sate = zeros(numel(sat),ne);
        for j=1:ne
            sate(:,j) = rSolE{j}.s;
        end
    end
    is = is + 1;
    ns(is) = numel(sat);    
end

%pressure
if ~isempty(strmatch('prf',variables))
    prf = rSol.pressure;
    if ne > 0
        prfe = zeros(numel(prf),ne);
        for j=1:ne
            prfe(:,j) = rSolE{j}.pressure;
        end
    end
    is = is + 1;
    ns(is) = numel(prf);
end

% structural parameters
if ~isempty(strmatch('stc',variables))
    fn = fieldnames(faults);
    dummy1 = [];
    for i = 1:length(fn)
        dummy2 = eval(['faults.' fn{i}]);
        dummy1 = [dummy1; reshape(dummy2,numel(dummy2),1)];
        is = is + 1;
        ns(is) = numel(dummy2);
    end
    stc = dummy1; clear dummy*
    if ne > 0 
        stce = zeros(numel(stc),ne);
        for j=1:ne
            dummy1 = [];
            for i = 1:length(fn)
                dummy2 = eval(['faultsE{j}.' fn{i}]);
                dummy1 = [dummy1; reshape(dummy2,numel(dummy2),1)];
            end
            stce(:,j) = dummy1; clear dummy*
        end
    end
end

xt = []; A = [];
for i = 1 : length(variables)
    xt = [xt; eval(variables{i})];
    if ne > 0
        A = [A; eval([variables{i} 'e'])];
    end
    dummy = ['clear ' variables{i} '*']; 
    eval(dummy);
end

nstat = ns; clear ns
