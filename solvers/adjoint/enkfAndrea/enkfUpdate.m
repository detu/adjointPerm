function [xam,Xap, Pa] = enkfUpdate(xfm, Xfp, Yt, H, R)
%% Inputs
% xfm: is the mean of the forecasted states
% Xfp: are the forecasted states
% Yt: is the actual measurement
% H: is the map between states and measurementes, Yt = H x

%% Outputs
% Xap: is the assimilated states Xap = Xfp + K [Y - H xf], 
% K = Pf H' [H Pf H' + Rk]^(-1)
% xam: is the assimilated mean of the states, considered as the best estimate

% Pa: is the error covariance matrix after the analysis step Pa = [I-K H] Pf

ny = size(H,1);
nx = size(H,2);
Ne = size(Xfp,2);

%identity matrix
Id = speye(nx,nx);

%IdMatr = speye(ny);
K2 = zeros(ny,ny);



% % reset random stream
% stream = RandStream.getGlobalStream
% reset(stream,1);

% generate perturbated measurements
Y = repmat(Yt,1,Ne);
Y = Y + normrnd(0,sqrt(R),ny,Ne);


% compute Pf
At = Xfp - Xfp* (ones(Ne,Ne)/Ne);
Pf = (At*At')/(Ne-1);


% compute kalman gain
K2 = inv(H *Pf *H' + R);
%[U,S,V] = svd(K2); % svd decomposition of  K2, K2 = U*S*V'
%Sdiag = di  ags(diag(1./S),0,ny,ny); %spdiags(diag(1./S),0,ny,ny)
K = Pf *H' * ( K2 + R ); % we can substitute R with the realized error covariance

% asimilate the states
Xap = Xfp + K * (Y - H *Xfp);

% compute the assimilated mean
xam  = xfm + K *(Yt -H*xfm);

% assimilated error covariance matrix 
Pa = (Id-K *H)* Pf;
