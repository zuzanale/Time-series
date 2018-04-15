function [vAlphasmooth,mV] = kf_smooth_st_NAIS(x,H,mT,c,d,sigmaEta,a0,P0) 

%% Kalman filter, matrix notation
vy = x; %v stands for vector, m for a matrix
T = size(vy,1);
vA = zeros(T,1);
mP = zeros(T,1);
vV = zeros(T,1);
mK = zeros(T,1);
vU = zeros(T,1);
mF = zeros(T,1);
mL = zeros(T,1);
mZ = 1;
mQ = sigmaEta;
mH = H;

mR = 1;
% initial values
a = a0;
p = P0;

%% filtering
for t = 1:T
    mP(t) = p;
    vA(t) = a;
    vU(t) = vA(t);
        
    % missing values threathment
       vV(t) = vy(t) - c - mZ*vA(t);
       mF(t) = mZ*mP(t)*mZ' + mH(t,t);
       mK(t) = mT*mP(t)*mZ'/(mF(t));
       mL(t) = mT - mK(t);
       a = mT*vA(t) + mK(t)*vV(t) + d;
       p = mT*mP(t)*mT' + mR*mQ*mR' - mK(t)*mF(t)*mK(t)';
end

%% loglikelihood evaluation due to (Q)MLE
l=  -(1/2)*T*log(2*pi) -(1/2)*sum(log(mF) +((vV.^2)./mF)); 
llik =mean(l);

%% Kalman smoothing
vR = zeros(T,1);
vAlphasmooth = zeros(T,1);
mN = zeros(T,1);
mV = zeros(T,1);

r = 0;
n = 0;
for t = T:-1:1
vR(t) = r;
mN(t) = n; 
    r = mZ/(mF(t))*vV(t) + mL(t)*r; %eq 4.38 pg.89
    vAlphasmooth(t) = vA(t) + mP(t)*r + d; %eq  4.39, first part

    n = mZ/(mF(t))*mZ' + mL(t)*n*mL(t)';%eq 4.42
    mV(t) = mP(t) - mP(t)*n*mP(t)' - c;%eq 4.40 or %eq 4.43 together
end
end