%% 2D 5-spot example with uncertain permeability.
%
% The model consists of a square horizontal reservoir section. Four
% producer sare positioned in the corners and an injector is positioned in
% the center. Permeability realizations are read from a file (see below).
%

%% model definition

    % Dimensions
    dims(1)=21; dims(2)=21; dims(3)=1; 
    physDims(1)=1000; physDims(2)=1000; physDims(3)=10;

    % Build grid
    grdecl.cartDims = reshape(dims, 1, []);
    [X, Y, Z]  = ndgrid(linspace(0, 1, dims(1) + 1), ...
                        linspace(0, 1, dims(2) + 1), ...
                        linspace(0, 1, dims(3) + 1));

    X = X .* physDims(1);
    Y = Y .* physDims(2);
    Z = Z .* physDims(3);

    % Make pilars
    lines          = zeros([prod(dims([1, 2]) + 1), 6]);
    lines(:,[1,4]) = reshape(X(:,:,[1,end]), [], 2);
    lines(:,[2,5]) = reshape(Y(:,:,[1,end]), [], 2);
    lines(:,[3,6]) = reshape(Z(:,:,[1,end]), [], 2);

    grdecl.COORD = reshape(lines', [], 1);

    % Assign z-coordinates
    ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
    z   = Z(ind(1), ind(2), ind(3));
    
    grdecl.ZCORN = z(:);

    % Assign active cells
    actnum = ones(dims);
    grdecl.ACTNUM = int32(actnum(:));

    % Process grid
    G = processGRDECL(grdecl);
    G = computeGeometry(G(1)); 
    clear grdecl

    % Rock properties
    %load([ROOTDIR 'examples\data\randomgauss21x21_spherical1000_conditioned.mat']);
    load('randomgauss21x21_spherical1000_conditioned.mat');
    
    for i=1:max(1,ne)
        perm(:,i) = K(:,i);
    end
    perm(:,ne+1) = K(:,end);
    K = perm; clear perm
    
%   Kmin=min(min(K)); Kmax=max(max(K));
%   K = (1 + (K-Kmin)/(Kmax-Kmin)*1000);
 
    K = 10 * exp(1 * K);
    
    phi = 0.25 * (K ./200).^0.1;
    
    K = convertFrom(K, milli*darcy);
    
    rock.perm = [K(:,end) K(:,end) K(:,end)/10];
    rock.poro = phi(:,end);
    
%   clear phi

    % Fluid properties [water, oil]
    fluid = initCoreyFluid('mu' , [   0.4,  0.9]*centi*poise  , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   4,   4]                 , ...
                           'sr' , [ 0.2, 0.2]                 , ...
                           'kwm', [   1,   1]);

    % Initial state
    p0 = 350*barsa();
    s0 = 0.2;
    
    % Well definitions and constraints
    W = [];
    I = 11;
    J = 11;
    q = 150*meter^3/day; %0.006*1000*meter^3/day;
    bhp = 350*barsa();
    nIW = 1:numel(I);
    for i = 1 : numel(I),
       W = verticalWell(W, G, rock, I(i), J(i), 1:dims(3), 'Type', 'rate', ...
                        'Val', q, 'Radius', 0.1, 'Comp_i', [1,0,0], ...
                        'name', ['I$_{', int2str(i), '}$']);           
    end

    I = [1  1 21 21];
    J = [1 21  1 21];
    bhp = 300*barsa();
    nPW = (1:numel(I))+max(nIW);
    for i = 1 : numel(I),
       W = verticalWell(W, G, rock, I(i), J(i), 1:dims(3), 'Type', 'bhp', ...
                        'Val', bhp, 'Radius', 0.1, ...
                        'name', ['P$_{', int2str(i), '}$']);
    end

    % gravity
    gravity off;                       
 
    