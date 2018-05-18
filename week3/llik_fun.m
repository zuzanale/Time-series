function [llik]=llik_fun(x,theta)            
%theta_ini = [sigma2_eta_ini; omega_ini; phi_ini]; %dimension 3x1
        a0=theta(2)/(1-theta(3));
        P0 = theta(1)/(1-theta(3)^2);
        H =(pi^2/2)*eye(size(x,1)); %variance of chi squared distribution
        sigmaEta= theta(1);
        c = -1.27;%mean 
        d = theta(2); %omega
        mT = theta(3); %phi     
    
%% Kalman filter derivation
  %% 1.initialization
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

    %% initial values
    a = a0;
    p = P0;

%% filtering
for t = 1:T
    mP(t) = p;
    vA(t) = a;
    vU(t) = vA(t);
        
       vV(t) = vy(t) - c - mZ*vA(t);
       mF(t) = mZ*mP(t)*mZ' + mH(t,t);
       mK(t) = mT*mP(t)*mZ'*inv(mF(t));
       mL(t) = mT - mK(t);
       a = mT*vA(t) + mK(t)*vV(t) + d;
       p = mT*mP(t)*mT' + mR*mQ*mR' - mK(t)*mF(t)*mK(t)';
end

%% loglikelihood evaluation due to (Q)MLE

l=  -(1/2)*T*log(2*pi) -(1/2)*sum(log(abs(mF)) +((vV.^2)./mF)); 
llik =mean(l);

end
      