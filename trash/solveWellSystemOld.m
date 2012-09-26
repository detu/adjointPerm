function [resSol, wellSol] =solveWellSystemOld(G, S, W)
%%
%% [resSol, wellSol] = solveWellSystem(Grid, System, Well)
%%
%% Solves System given fields BI,C,D RHS.
%%
%%   resSol.pressure   : pressure in active cells
%%   resSol.localFlux  : out-flux over local faces
%%   resSol.globalFlux : flux over global indexed faces  corr to G.faces.neoghbors
%%   wellSol.flux
%%   wellSol.pressure


%% Build wellsystem
numWells = length(W);
BIW = []; CW = []; DW = [];
fW = []; hW = [];
inx = struct;
for k = 1:numWells
  w = W(k);
  sb1 = size(BIW); sb2 = w.S.sizeB;
  BIW = [BIW            sparse(sb1(1), sb2(2))
        sparse(sb2(1), sb1(2))  w.S.BI];
  inx(k).flux = (1:length(w.S.RHS.f)) + length(fW);
  fW = [fW; w.S.RHS.f];

  CW = [CW
        w.S.C];

  sd1 = size(DW); sd2 = w.S.sizeD;
  DW = [DW            sparse(sd1(1), sd2(2))
        sparse(sd2(1), sd1(2))  w.S.D];
  inx(k).pres = (1:length(w.S.RHS.h)) + length(hW);
  hW = [hW; w.S.RHS.h];
end

%% Build full system
BI = [S.BI           sparse( size(S.BI,1), size(BIW,2))
      sparse( size(BIW,1), size(S.BI,2))  BIW   ];

C = [S.C
     CW];

D = [S.D       sparse( size(S.D,1), size(DW,2) )
     sparse( size(DW,1), size(S.D,2))  DW ];

f = [S.RHS.f; fW];
g = S.RHS.g;
h = [S.RHS.h; hW];

%% Solution

l=diag(C'*BI*C);
LI=spdiags(1./l,0,length(g),length(g));

disp('Computing matrix products etc ...'); tic
M = C'*BI*D;
Sym = D'*BI*D;

fHat = g - C'*(BI*f);
gHat = h - D'*(BI*f);
toc

disp('Solving linear system:'); tic
%if length(gHat)<100000
  lam=(Sym-M'*LI*M)\(M'*(LI*fHat)-gHat);
%else
%  disp('amg:')
%  lam=amg(Sym-M'*LI*M, M'*(LI*fHat)-gHat,zeros(size(gHat)), ...
%          10^(-10),200,10,10,5);
%end
toc

pres = LI*fHat+LI*(M*lam);
flux = BI*(f-D*lam+C*pres);


solution.pressure = pres;
solution.lambda = lam;

%% reservoir solution
resSol.pressure = pres;
resSol.lambda = lam(1:S.sizeD(2));
resSol.localFlux = flux(1:S.sizeB(1));
resSol.globalFlux = cellFlux2faceFlux(G,resSol.localFlux);

%% well solution
pW = lam(S.sizeD(2)+1:end);
fW = flux(S.sizeB(1)+1:end);

wellSol = struct;
for k=1:numWells
  wellSol(k).flux = -fW(inx(k).flux);
  wellSol(k).pressure = pW(inx(k).pres);
end

return
