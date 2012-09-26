function partitionNU = partitionNonuniform(lower_bound,upper_bound,G,rock,rSol,varargin)
%
% SYNOPSIS:
%   partitionNU = partitionNonuniform(lower_bound, upper_bound,G,rock,rSol)
%   partitionNU = partitionNonuniform(lower_bound, upper_bound,G,rock,rSol,...
%                                     'Wells',W,'wSol',wSol)
%   partitionNU = partitionNonuniform(lower_bound, upper_bound,G,rock,rSol,'src',src)
%   
%
%
% DESCRIPTION:
%   Partitions the grid into a non-uniform grid partitioning according to
%   the non-uniform grid coarsening algorithm in
%     J. E. Aarnes, V. L. Hauge, and Y. Efendiev:
%     Coarsening of three-dimensional structured and unstructured grids for subsurface flow.
%     Advances in Water Resources, Vol. 30, Issue 11, November 2007, pp. 2177-2193.
%     DOI: 10.1016/j.advwatres.2007.04.007
%
%    OBS: Subfunksjonene her er ikke pyntet ferdig. Fortsatt under
%    utvikling...
%
% PARAMETERS:
%   upper_bound
%   lower_bound
%   G  - Grid data structure.
%   rock - Rock structure.
%   rSol - Reservoir solution structure. Uses the field faceFlux.
%   
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells  - Well structure as defined by function 'addWell'.  May be empty
%            (i.e., W = []) which is interpreted as a model without any
%            wells.
%
%   wSol   - Well solution data structure.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = []) which is
%            interpreted as a reservoir model without explicit sources.
% 
%   square - Default true: Parameters that determines if the coarse blocks
%            have square shape or diamond shape (original shape, use false).
%
%   plotif - Parameter that determines if we plot the steps in the
%            algorithm. Default false.
%
% EXAMPLE:
%        [G, W, rock] = SPE10_setup(2);
%        rock.poro = max(rock.poro, 1e-3);
%
%        fluid  = initSimpleFluid('rho', [1000.0, 800.0], 'mu', [1.0, 1.0]);
%
%        rSol = initResSol(G, 0.0, 0.0);
%        wSol=initWellSol(W,0.0);
%        S = computeMimeticIP(G, rock);
%        [rSol,wSol]= solveIncompFlow(rSol, wSol, G, S, fluid,'wells', W);
%
%        lower_bound=10; upper_bound=50;
%        partitionNU=partitionNonuniform(lower_bound,upper_bound,G,rock,rSol, 'Wells',W,'wSol',wSol, 'plotif', true);
%        CG_NUC=generateCoarseGrid(G,partitionNU);
%
%        %% call on transport solver:
%        rSol  = implicitTransportCG(rSol, wSol, G, CG_NUC, dT, rock, fluid, 'wells', W);
%
%
% RETURNS:
%   partitionNU - A partition vector of non-uniformly coarsened grid.
%
% SEE ALSO:

%{
#COPYRIGHT#
%}

% $Date: 2009-10-26 09:19:44 +0100 (ma, 26 okt 2009) $
% $Revision: 3067 $

opt = struct('Wells',[],'wSol',[],'src', [], 'square',true, 'plotif',false);
opt = merge_options(opt, varargin{:});

plotif=opt.plotif;

%% Compute a representative flux value for the cell center.
cellFlux = computeCellFlux(G, rock, rSol, 'Wells', opt.Wells, 'wSol', opt.wSol,'src',opt.src); 
if(plotif)
   figure; plotCellData(G, log10(cellFlux')); axis equal tight off;
   title('Log of cellFlux');
end


%% The four steps in the algorithm
%% Generate initial partitioning of the velocities into 10 bins
partitionNU=partitionInitialNUC(G, cellFlux);
if(plotif)
   figure;
   subplot(2,2,1); outlineCoarseGrid(G, partitionNU); hold on; outlineDomain(G);
   axis equal tight off;
   title(sprintf('Initial partitioning of cells: %d',max(partitionNU)));
end

%% Merge too small blocks with neighboring blocks
partitionNU=mergeBlocks(G, partitionNU,cellFlux, lower_bound,rock);
if(plotif)
   subplot(2,2,2); outlineCoarseGrid(G, partitionNU); hold on; outlineDomain(G);
   axis equal tight off;
   title(sprintf('After merging step: %d',max(partitionNU)));
end

%% Refine too large blocks
partitionNU=refineBlocks(G,partitionNU,cellFlux,upper_bound,rock,'square',opt.square);
if(plotif)
   subplot(2,2,3); outlineCoarseGrid(G, partitionNU); hold on; outlineDomain(G);
   axis equal tight off;
   title(sprintf('After refining step: %d',max(partitionNU)));
end

%% Merge too small blocks with neighboring blocks
partitionNU=mergeBlocks(G,partitionNU,cellFlux, lower_bound, rock);
if(plotif)
   subplot(2,2,4); outlineCoarseGrid(G, partitionNU); hold on; outlineDomain(G);
   axis equal tight off;
   title(sprintf('After mergeing step: %d',max(partitionNU)));
end


end

function cellFlux = computeCellFlux(G, rock, rSol, varargin)
% SYNOPSIS:
%   cellFlux = computeCellFlux(G, rock, rSol)
%   cellFlux = computeCellFlux(G, rock, rSol, 'src', src)
%   cellFlux = computeCellFlux(G, rock, rSol, 'Wells', W, 'wSol', wSol)
%
%
% DESCRIPTION:
%   Computes a representative scalar value of the flux in each cell in the
%   grid.
%
% PARAMETERS:
%   G  - Grid data structure.
%   rock - Rock structure
%   rSol - Reservoir solution structure. Uses the field faceFlux.
%
%
% RETURNS:
%   cellFlux - A vector with a scalar value of the flux in each grid cell.
%
% SEE ALSO:
%   partitionNonuniform, mergeBlocks, refineBlocks

opt = struct('Wells',[],'wSol',[],'src', []);
opt = merge_options(opt, varargin{:});

pv=poreVolume(G,rock);

V=zeros(1,G.cells.num);

%% Compute a representative value of fluxes in each cell
% For all cells we find the faces belonging to this cell. The flux on these
% faces is multiplied with the distance from the cell centroid to the face
% centroid and summed. The total sum is divided by the pore volume of the
% cell.
for i=1:G.cells.num,
   ei = [G.cells.facePos(i):1:G.cells.facePos(i+1)-1];
   ef = G.cellFaces(ei, 1);
   v_i=0;
   for j=1:length(ef),
      % For all faces, sum the faceflux times distances edge to cell centroid.
      v_i=v_i+rSol.faceFlux(ef(j))*(G.faces.centroids(ef(j),:)-G.cells.centroids(i,:));
   end;
   V(i)=norm(v_i/pv(i)); 
   %V(i)=norm(v_i/G.cells.volumes(i));
end
%{
%  Equivalent:
cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2).';
C      = G.faces.centroids(G.cellFaces(:,1), :) - ...
         G.cells.centroids(cellNo,           :);
val    = abs(rSol.faceFlux(G.cellFaces(:,1))) .* sqrt(sum(C .* C, 2));
V      = accumarray(cellNo, val, [G.cells.num, 1]) ./ G.cells.volumes;
%}

%% Source/sink cells and well cells need special consideration
%  For source/sink cells we use the maximum value of the face fluxes and
%  divide by the area of this face. 
%  For well cells we use the flux from the well solution struct for this
%  cell divided by the sum of areas of the faces belonging to this cell.
%  In both cases we divide by the porosity of the cell.
if ~isempty(opt.src)
   sourceCells=opt.src.cell; % both sources and sinks
   
   for i=1:length(sourceCells),
      cell=sourceCells(i);
      ei = [G.cells.facePos(cell):1:G.cells.facePos(cell+1)-1]; % For raskere kode
      ef = G.cellFaces(ei, 1);
      [val c]=max(abs(rSol.faceFlux(ef)));
      V(cell)=val/G.faces.areas(ef(c))/rock.poro(cell);
   end
   %{
   % Equivalent
   function r = maxloc(x),
      [m, i] = max(x);
      r = {[m, i-1]};
   end
   c  = opt.src.cell;
   nc = numel(c);
   [s, e] = deal(G.cells.facePos(c), G.cells.facePos(c + 1) - 1);
   nf = diff([s, e], [], 2) + 1;
   m  = accumarray(rldecode((1 : nc) .', nf),                       ...
                   abs(rSol.faceFlux(G.cellFaces(mcolon(s,e), 1))), ...
                   [nc, 1], @maxloc);
   m  = vertcat(m{:});
   V(c) = m(:,1) ./ G.faces.areas(G.cellFaces(s + m(:,2), 1));
   %}
end

% Wells cells
% wSol(celle)/SUM(areal av facene)
% når well celle midt i grid skal være lik max(rSol(ef))/areal(ef)
if ~isempty(opt.Wells)
    for i=1:length(opt.Wells)
        for j=1:length(opt.Wells(i).cells)
            wellCell=opt.Wells(i).cells(j);
            ei = [G.cells.facePos(wellCell):1:G.cells.facePos(wellCell+1)-1]; % For raskere kode
            ef = G.cellFaces(ei, 1);
            area=sum(G.faces.areas(ef));
            V(wellCell)=abs(opt.wSol(i).flux(j))/area/rock.poro(wellCell);
        end
    end
    %{
    % Equivalent:
    c  = vertcat(opt.Wells.cells);
    nc = numel(c);
    [s, e] = deal(G.cells.facePos(c), G.cells.facePos(c + 1) - 1);
    nf = diff([s, e], [], 2) + 1;
    V(c) = abs(vertcat(opt.Wells.flux)) ./ ...
           accumarray(rldecode((1 : nc).', nf), ...
                      G.faces.areas(G.cellFaces(mcolon(s,e), 1));
    %}
end

cellFlux = V;
%{
% Ad hoc hack... for trekantgriddet... unngå 0'er
c0=find(cellFlux==0);
c=find(cellFlux);
[val idx]=min(cellFlux(c));
cellFlux(c0)=val;
%}
end


function partitionNU = partitionInitialNUC(G, cellFlux)
%
% SYNOPSIS:
%   partition = partitionInitialNUC(G, cellFlux)
%
% DESCRIPTION:
%   Generates the initial partitioning of a non-uniform coarse grid.
%   The grid is partitioned into ten velocity groups according to the
%    logarithm of the velocity in the grid blocks. 
%
% PARAMETERS:
%   G  - Grid data structure.
%   cellFlux - vector of a representative flux for each grid cell.
%
% RETURNS:
%   partitionNU - a partition vector 
%
% SEE ALSO:
%   computeCellFlux, mergeBlocks, refineBlocks

partition = zeros(G.cells.num,1);

Li=10; Lg=log10(cellFlux);
Mg=max(Lg)+1e-5;
mg=min(Lg)-1e-5; dl=abs(Mg-mg)/Li;
%
%% Initial non-uniform grid based on flow information
% In ten steps we divide the grid cells into bins with similar logarithm
% of velocity
for i=1:Li
    c=find(Lg>=mg+(i-1)*dl & Lg<mg+i*dl);
    partition(c)=i;
end

partitionNU = compressPartition(partition);
partitionNU = processPartition(G, partitionNU, 'Reconnect', false);
partitionNU = compressPartition(partitionNU);

end


function partition = mergeBlocks(G, partition,cellFlux, lower_bound, rock, varargin)
%
% SYNOPSIS:
%   partition = mergeBlocks(G,partition,cellFlux,lower_bound,rock);
%   partition = mergeBlocks(G,partition,cellFLux,lower_bound, rock, ...
%                           'reorder', reorderpartition);
%
%
% DESCRIPTION:
%     Merges blocks with volume less than "lower_bound"/N fraction
%     of the total resevoir volume with a neighboring block subject to
%     flow of similar magnitude.
%     OBS: under utvikling...
%
% PARAMETERS:
%   G  - Grid data structure.
%   partition - an initial partitioning of the grid according to the
%               logarithm of the velocity. (The output of partitionNonuniform).
%   cellFlux - A vector of scalar values of the flux/velocity in each grid
%              cell.
%   lower_bound - the coarsening parameters giving the lower bound of the
%                 coarse grid blocks.
%   rock        - Rock structure
%
%
% RETURNS:
%   partition - partition vector of the grid.
%
% SEE ALSO:
%   computeCellFlux, partitionNonuniform, refineBlocks

opt = struct('reorder',[]);
opt = merge_options(opt, varargin{:});
% Denne reorder opsjonen skal vel kanskje bort etterhvert...
% Ikke ferdig i alle fall... derfor er den ikke med i ytterste kall...

pv=poreVolume(G,rock);

% if combining with reorder...
tagged_grid=opt.reorder;

%fluxVol=log10(cellFlux)'.*G.cells.volumes;
fluxVol=log10(cellFlux)'.*pv;
%lower_bound_volume = lower_bound*sum(G.cells.volumes)/G.cells.num;
lower_bound_volume = lower_bound*sum(pv)/G.cells.num;

%blockVolumes=accumarray(partition,G.cells.volumes);
blockVolumes=accumarray(partition,pv);
blockFluxVol=accumarray(partition,fluxVol);

da=blockFluxVol./blockVolumes; % Finne bedre navn

skal_slettes=zeros(length(partition),1);

%%  Merging algorithm
% As long as we have blocks with too small volumes, we attempt to find a
% neighbor with similar flow magnitude and merge the block to this block.

while(min(blockVolumes) < lower_bound_volume)
   
   % Fra connections in genreateCoarseGrid.m
   nb=max(partition);
   p1=[0; partition];
   pN=p1(G.faces.neighbors+1);
   
   if(isfield(G.faces,'tag')),
      pN=pN(find(~G.faces.tag),:); % exlude faces that are faults
   end
   pN=pN(all(pN>0,2),:);
   
   %pN = pN(pN(:,1) ~= pN(:,2), :);% Exclude block-internal connections.
   
   internal_idx= find(pN(:,1) ~= pN(:,2));
   if(isfield(G.faces,'tag'))
      internal_notfault= setxor(internal_idx, find(G.faces.tag));
      pN=pN(internal_notfault,:);
   else
      pN=pN(internal_idx,:);
   end
   
   % Build (symmetric) connectivity matrix.
   EE = sparse([pN(:,1); pN(:,2); (1 : nb).'], ...
      [pN(:,2); pN(:,1); (1 : nb).'], 1);   %OBS denne kalles mange ganger og tar tid... Alternativ?
   
   % Find the block with smallest volume
   %[val block]=min(blockVolumes);
   
   % Find blocks with too small volumes
   block_array=find(blockVolumes < lower_bound_volume);
   
   for i=1:length(block_array),
      block=block_array(i);
      cellsInBlock=find(partition==block); %CG.cells.subCells(:,block);
      
      ei=EE(block,:); % find faces belonging to this block
      ei(block)=0;     % fjerner seg selv...
      naboer=find(ei); % ci1
      
      %   [cpos f]=findBlockFaces(G,partition);
      %   blockFaces = f(cpos(block) : cpos(block+1) - 1);
      
      noMoreToDo=0;
      if(isempty(naboer))
         % Hvis ikke har noen naboer; typisk kan skje når har faults som hindrer
         % merging
         skal_slettes(cellsInBlock)=1;
         noMoreToDo=1;
      end;
      
      if(~noMoreToDo)
         di1=abs(da(naboer)-da(block));
         ci=find(abs(di1-min(di1))<1e-12);
         %[val i]=min(abs(da(naboer)-da(block)));
         [mergeToBlock idx]=min(naboer(ci));
         
         cellsInNewBlock=find(partition==mergeToBlock); %CG.cells.subCells(:,mergeToBlock);
         found=1;
         % Make sure we merge block with a block with same reorder tag.
         if(~isempty(tagged_grid))
            tag=tagged_grid(cellsInBlock); tag=tag(1);
            tag2=tagged_grid(cellsInNewBlock); tag2=tag2(1);
            
            found=0; teller=1;
            if(tag==tag2) found=1; end
            
            naboer_2=naboer;
            while(~found & teller<nnz(naboer_2))
               naboer_2=setxor(naboer_2,mergeToBlock); % Er setexor dyr?
               [mergeToBlock idx]=min(naboer_2); % finne ny nabo
               
               cellsInNewBlock=find(partition==mergeToBlock); %CG.cells.subCells(:,mergeToBlock);
               tag2=tagged_grid(cellsInNewBlock); tag2=tag2(1);
               
               if(tag2==tag) found=1; end
               teller=teller+1;
            end
         end
         if(found)
            partition(cellsInBlock)=mergeToBlock; % "Sprang" i partisjonsvektoren nå...
            
            %Oppdatere EE etter Jørgs kode:
            EE(mergeToBlock, naboer)=EE(mergeToBlock, naboer) + EE(block, naboer);
            EE(naboer, mergeToBlock)=EE(naboer, mergeToBlock) + EE(naboer, block);
            EE(:,block)=0;
         else
            skal_slettes(cellsInBlock)=1;
         end
      end
      
   end
   partition=compressPartition(partition);
   %partition=processPartition(G,partition,'Reconnect', false);
   
   % Finne nye verdier
   %V=G.cells.volumes;
   V=pv;
   V(find(skal_slettes))=lower_bound_volume*1000;
   blockVolumes=accumarray(partition,V);
   blockFluxVol=accumarray(partition,fluxVol);
   da=blockFluxVol./blockVolumes;
   
end

end


function partition = refineBlocks(G,partition, cellFlux,upper_bound, rock, varargin)
%
% SYNOPSIS:
%   partition = refineBlocks(G,partition,cellFlux,upper_bound,rock);
%   partition = refineBlocks(G,partition,cellFLux,upper_bound,rock, ...
%                           'square', true);
%
%
% DESCRIPTION:
%     Refine to large blocks/blocks with too large flux.
%
%     OBS: under utvikling...
%
% PARAMETERS:
%   G  - Grid data structure.
%   partition - an initial partitioning of the grid according to the
%               logarithm of the velocity. (The output of partitionNonuniform).
%   cellFlux - A vector of scalar values of the flux/velocity in each grid
%              cell.
%   upper_bound - the coarsening parameters giving the upper bound on the flow
%                 through each coarse block.of the
%   rock - rock structure
%                 
%
% RETURNS:
%   partition - partition vector of the grid.
%
% SEE ALSO:
%   computeCellFlux, partitionNonuniform, mergeBlocks

%
opt = struct('square',false);  % true or false
opt = merge_options(opt, varargin{:});

pv=poreVolume(G,rock);

%dv=(log10(cellFlux)-min(log10(cellFlux))+1)'.*G.cells.volumes;
dv=(log10(cellFlux)-min(log10(cellFlux))+1)'.*pv;
Va=upper_bound*sum(dv)/G.cells.num;

% Naboskap
Neighbors = G.faces.neighbors;
N=G.cells.num;

ni=find(Neighbors(:,1)); Neighbors=Neighbors(ni,:);  % remove zeros...
ni=find(Neighbors(:,2)); Neighbors=Neighbors(ni,:);
e1=double(Neighbors(:,1))'; e2=double(Neighbors(:,2))';
NeighborMatrix=sparse([1:N,e1,e2],[1:N,e2,e1],1);

numBlocks=max(partition); % counter for making new blocks

%% Refining algorithm
%  Go through each coarse block in the partition and refine if necessary.
for i=1:max(partition)
    cellsInBlock=find(partition==i);
    numCellsInBlock=length(cellsInBlock);
    Ci=G.cells.centroids(cellsInBlock,:)';
    Ci=Ci-repmat(Ci(:,1),1,numCellsInBlock);
    dc=sum(abs(Ci));  % distance vector
    cellCount=0; % counting cells in new block
    
    while cellCount<numCellsInBlock
        [maxVal,cellIdx]=max(dc); cellIdx=min(cellIdx);
        [idx4cellsInBlock,jj,cells]=find(cellsInBlock); %l=length(fr);
        e=zeros(length(idx4cellsInBlock),1);
        e(find(idx4cellsInBlock==cellIdx))=1;
        Ve=0; Veo=-1;
        GA=NeighborMatrix(:,cells);I=GA(cells,:);
        
        while Ve>Veo & Ve<Va
            e=I*e;  %find neighboring cells
            
            % making square blocks, not diamond shaped. (add diagonal
            % neighbors)
            if(opt.square)
                % ----- 
                % For å lage firkantstruktur, ikke diamant:
                e=spones(e);
                e_neigh=I*e;  % finner naboenes naboer
                e_neigh=spones(e_neigh); %e_neigh=e_neigh>=1;
                e_neigh=e_neigh-e;%trekker ut de "nye" naboene
                fe_neigh=find(e_neigh); % angir indeksene til nye naboer
                % finner hvilke av de nye naboene som er nabo til 2 av cellene i e:
                for j=1:size(fe_neigh)
                    no_neigh=I(fe_neigh(j),:)*e; % Finner hvor mange felles naboer
                    % no_neigh har med cellene i e
                    if(no_neigh==1)
                        e_neigh(fe_neigh(j))=0;  % Hvis bare 1 felles, slettes denne cellen
                        % hvis det er mer enn 1 felles beholdes den
                    end
                end
                e=e+e_neigh;  % Legger til nye "naboer" i tillegg til de
                % opprinnelige i e.
                % -----
            end
            
            Veo=Ve+1e-15;
            neighboringCells=find(e); ce=cells(neighboringCells);
            Ve=sum(dv(ce));
        end
        numBlocks=numBlocks+1;
        partition(ce)=numBlocks;  % Add new block in the partition vector
        cellsInNewBlock=idx4cellsInBlock(neighboringCells);
        dc(cellsInNewBlock)=0;  % Update the distance vector
        cellsInBlock(cellsInNewBlock)=0; % Remove these cells as they no longer are counted in the original block
        cellCount=cellCount+1;
    end
end

% Update the partition vector
%partition=processPartition(G,partition,'Reconnect', false);
partition=compressPartition(partition);

end



function varargout = outlineDomain(G)
%Impose outline of reservoir domain/outer edges of grid on existing grid plot.
%
% SYNOPSIS:
%       outlineDomain(G)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.
%
% EXAMPLE:
%   G  = cartGrid([8, 8, 2]);
%   p  = partitionUI(G, [2, 2, 1]);
%   % plot fine grid:
%   plotGrid(G, 'faceColor', 'none'); view(3);
%   % outline coarse grid on fine grid:
%   outlineCoarseGrid(G, p);
%   outlineDomain(G);
%
% SEE ALSO:
%   plotFaces, plotGrid, patch.

nx=G.cartDims(1); ny=G.cartDims(2); nz=G.cartDims(3);
% 2D
if(G.cartDims(3)==1)
    vertices=[1,nx+1,(nx+1)*(ny+1), (nx+1)*ny+1, 1];
    pts=G.nodes.coords(vertices,:);
    plot(pts(:,1), pts(:,2),'k');
    % Viser ikke alle linjene i figuren. Bug i matlab...
end

end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function [fp, f] = findBlockFaces(G, p)
   p1 = [0; p];
   pN = p1(G.faces.neighbors + 1);
   fn = (1 : size(G.cellFaces, 1)) .';
   i  = pN(:,1) ~= pN(:,2);
   i  = i(G.cellFaces(:,1));
   cn = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';

   f  = sortrows([p(cn(i)), fn(i)]);
   fp = cumsum([1; accumarray(f(2:end,1), 1, [max(p), 1])]);
   f  = f(2:end,2);
end

