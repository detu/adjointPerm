%% 2D example with structural uncertainty.
%
% The model consists of a horizontal reservoir cross section between
% vertical injector and producer wells. A sloping fault with non-zero throw
% is positioned somewhere between the 2 wells.
%

%% ONLY FOR TESTING
%   ne = 1; j = 3; S = repmat([0 0.5 0]',1,10) + repmat([0.1 0.1 0.2]',1,10) .* randn(3,10);

%% model definition

    % Dimensions
    dims(1)=100; dims(2)=1; dims(3)=3;
    physDims(1)=1000; physDims(2)=20; physDims(3)=100;

    % Build grid
    grdecl.cartDims = reshape(dims, 1, []);
    [X, Y, Z]  = ndgrid(linspace(0, 1, dims(1) + 1), ...
                        linspace(0, 1, dims(2) + 1), ...
                        linspace(0, 1, dims(3) + 1));

    X = X .* physDims(1);
    Y = Y .* physDims(2);
    Z = Z .* physDims(3); 
    
    % fault slope
    if isempty(S)
        a1 = 0.1 * randn(1,ne);
    else
        a1 = S(1,:);
    end
    if j == 0
        a = -0.1;
    else
        a = a1(1,j);
    end
    
    % fault position
    if isempty(S)
        b1 = 0.5 + 0.1*randn(1,ne);
    else
        b1 = S(2,:);
    end
    if j == 0
        b = 0.35;
    else
        b = b1(1,j);
    end

    if b < 0.025 - a, b = 0.025 - a + b; end
    if b > 0.975 + a, b = 2 * (0.975 + a) - b; end
    
    a = a * physDims(1);
    b = b * physDims(1);  
    
    jj=linspace(0, 1*physDims(1), dims(1) + 1); jb=find(jj<=b,1,'last');
    
    dummy1 = find(X<=b); dummy2=find(X>b);  
    X(dummy1)=X(dummy1)-a*(X(dummy1)/b).*(Z(dummy1)-0.5*physDims(3))/physDims(3);
    X(dummy2)=X(dummy2)-a*((1*physDims(1)-X(dummy2))/(1*physDims(1)-b)).*(Z(dummy2)-0.5*physDims(3))/physDims(3);
    clear dummy*

    %% Make pilars
    lines          = zeros([prod(dims([1, 2]) + 1), 6]);
    lines(:,[1,4]) = reshape(X(:,:,[1,end]), [], 2);
    lines(:,[2,5]) = reshape(Y(:,:,[1,end]), [], 2);
    lines(:,[3,6]) = reshape(Z(:,:,[1,end]), [], 2);

    grdecl.COORD = reshape(lines', [], 1);

    %% Assign z-coordinates
    % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
    ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
    z   = Z(ind(1), ind(2), ind(3));

    % fault throw
    if isempty(S)
        c1 = 0.15*randn(1,ne);
    else
        c1 = S(3,:);
    end
    if j == 0
        c = -0.15;
    else
        c = c1(1,j);
    end
    z(2*(jb-1)+1:end,:,:) = z(2*(jb-1)+1:end,:,:) + c * physDims(3);  
    
    grdecl.ZCORN = z(:);

    %% Assign active cells
    actnum = ones(dims);
    grdecl.ACTNUM = int32(actnum(:));

    %% Process grid
    G = processGRDECL(grdecl);
    G = computeGeometry(G(1)); 
    clear grdecl

    S(1,:) = a1; % / physDims(1);
    S(2,:) = b1; % / physDims(1);
    S(3,:) = c1;
    
    %% Rock properties
    K = 20 *ones(dims(1)*dims(2)*dims(3),ne+1);
    for i=1:ne+1
        dummy = reshape(K(:,i),dims(1)*dims(2),dims(3));
        dummy(:,1) = 2000; dummy(:,3) = 2000;
        K(:,i) = reshape(dummy,dims(1)*dims(2)*dims(3),1);
    end
    
    phi = 0.25 * (K ./200).^0.1;
    K = convertFrom(K,milli*darcy);

    rock.perm = [K(:,end) K(:,end) K(:,end)/10]; 
    rock.poro = phi(:,end);

    %% Fluid properties
    fluid = initCoreyFluid('mu' , [   0.4,  0.9]*centi*poise, ...
                         'rho', [1014, 859]*kilogram/meter^3, ...
                         'n'  , [   4,   4]                 , ...
                         'sr' , [ 0.2, 0.2]                 , ...
                         'kwm', [   1,   1]);

    %% Initial state
    p0 = 350*barsa();
    s0 = 0.2;
    
    bhp = 300*barsa();
    q = 1*0.003*1000*meter^3/day;

    %% Wells
    
    % Set vertical injectors
    I = [1];
    J = [1];
    R = q;
    bhp = 350.00*barsa(); P = bhp;
    nIW = 1:numel(I); W = [];
    for i = 1 : numel(I),
       W = verticalWell(W, G, rock, I(i), J(i), 1:dims(3), 'Type', 'bhp', ...
                        'Val', P(i), 'Radius', 0.1, 'Comp_i', [1,0,0], ...
                        'name', ['I$_{', int2str(i), '}$']);                
    end

    % Set vertical producers
    I = [dims(1)];
    J = [1];
    bhp = 340*barsa(); P = bhp;
    nPW = (1:numel(I))+max(nIW);
    for i = 1 : numel(I),
       W = verticalWell(W, G, rock, I(i), J(i), 1:dims(3), 'Type', 'bhp', ...
                        'Val', P(i), 'Radius', 0.1, ...
                        'name', ['P$_{', int2str(i), '}$']);
    end
    
    % ONLY FOR TESTING
%     figure(j+1)
%     axis equal
%     plotCellData(G,convertTo(rock.perm(:,1),milli*darcy),'EdgeColor','k','FaceAlpha', 1.0);
%     plotWell(G, W, 'height', 75, 'color', 'k');
%     axis off, view(15,60), h=colorbar('horiz');
%     caxis([0 2000]);
%     zoom(2.5), title('permeability [mD]');
%     axis tight off, view([-40, 30]);
    
  