function [llik]=llik_fun_IS(x,theta,eta_plus,u_plus)     
%debuggin
   %theta(1)=sigma2_u_ini;
    %theta(1)=sigma2_eta_ini ;
    %theta(2)=omega_ini; 
    %theta(3)=phi_ini ;
    %x=y;
   % eta_plus=eta_plus;
   % u_plus=u_plus;
    
%theta_ini = [sigma2_eta_ini; omega_ini; phi_ini]; %dimension 3x1
        a0=theta(2)/(1-theta(3));
        P0 = theta(1)/(1-theta(3)^2);
        %H =theta(1)*eye(size(x,1)); %this one might not be necessary to compute
        sigmaEta= theta(1);
        c = -1.27;%mean 
        d = theta(2); %omega
        mT = theta(3); %phi  
        
        %for conveniece, since somewhere we use sigma and omega and phi not
        %d or c or theta(x)
        sigma2_eta=theta(1);
        omega=theta(2);
        phi=theta(3);
%% Kalman filter derivation
  %% 1.initialization
    vy = x; %v stands for vector, m for a matrix
    T = size(vy,1);
    I = ones(T,1);
    g =  zeros(T,1);
    n=1;
    
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
    
%% MODE estimation
if g(n,1)==0 
    A = diag( (2 .*exp(g))./vy.^2 );
    vA = (2 .*exp(g))./vy .^2;
    z = g - (1/2).*vA+ I ;  
    Yn = z;
    H = A;    
[~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
g = modsmooth;
 n=n+1;
else
    while sqrt(g(n,1)-g(n-1,1)^2)>0.00001 % for cycle is more suitable than while cycle; after 15 iterations convergence is reached surely  
        %while cycle
        A = diag( (2 .*exp(g))./vy.^2 );
        vA = (2 .*exp(g))./vy .^2;
        z = g - (1/2).*vA+ I ;  
        Yn = z;
        H = A;    

    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
    g = modsmooth;
    n=n+1;
    end
end

%% Step 1. Construct approximating g
  A = diag( (2 .*exp(modsmooth))./y.^2 );%matrix
  vA = (2 .*exp(modsmooth))./y .^2;%vector
  z = g - (1/2).*vA+ I ;  
  %based on chapter 14 in DK book about y_t close to 0
  for t=1:T
      if z(t,1)<0.0000001
           z(t,1)=0.001;
      else
          z(t,1) = z(t,1);
      end
  end
  y_star = z;
  H = A;    
  
  %% Step 2. Evaluate log g(Yn)
llik_g =llik_fun(y_star,theta_hat); %evaluate log[g(Yn)]

%% Step 3 Construct log w
%% 3a) Draw theta^(i)
N = 100;
y_plus = zeros(T,1);
theta_plus= zeros(T+1,1);
theta_tilde = zeros(T,N);
p_Y_theta = zeros(1,N);
g_Y_theta = zeros(1,N);
m = zeros(1,N);
w = zeros(1,N);
a0 = omega / (1 - phi);
P0 = sigma2_eta / (1 - phi^2);
for i = 1:N 

% KFS to obtain conditional theta_hat
[~,con_theta_hat] =kf_smooth(y_star,A,phi,c,omega,sigma2_eta,a0,P0);

% KFS to obtain uncoditional simulation theta^+
    alpha_plus(1,1)=a0;%initialize alpha^+
    theta_plus(1,1)=alpha_plus(1,1);%initialize theta+
    for t = 1:T
        %u_plus = sqrt(vA(t))*randn(1,1);  % generate a vector discrepances
        %eta_plus = sqrt(sigma2_eta)*randn(1,1); %generate a vector discrepances 
        alpha_plus(t+1,1)=omega + phi * alpha_plus(t,1)+sqrt(sigma2_eta)*eta_plus;
        theta_plus(t+1,1)=alpha_plus(t+1,1);% i dont think this one is necessary
        y_plus(t,1)=  theta_plus(t,1)+ c +u_plus;
    end

% the smoothed mean from the simulate linear gaussian model
[~,con_theta_hat_plus] =kf_smooth(y_plus,A,phi,c,omega,sigma2_eta,a0,P0);

% uconditional simulation smoothing
theta_tilde(:,i) = con_theta_hat + theta_plus(1:end-1,1) -con_theta_hat_plus; %since the last value of vh_P is for t+1

%% Step 3b Evaluate p(Yn|theta) and g(Yn|theta)
%compute p(Yn|theta)
lp_Y_theta=-0.5.*(log(2*pi*exp(0.5*omega/(1-phi))^2) + theta_tilde(:,i) + z.^2.*exp(-theta_tilde(:,i)));
p_Y_theta(1,i)=exp(sum(lp_Y_theta,1));%numerically more stable version
%p_Y_theta(1,i)=prod(exp(lp_Y_theta),1);

%compute g(Yn|theta)
lg_Y_theta=-0.5.*(log(2*pi*vA)) - (1./(2*vA)).*(y_star-c-theta_tilde(:,i)).^2;%not sure if c supposted to be here
%g_Ytheta=1/(sqrt(2*pi*vA)).*exp(-((y_star-theta_tilde(:,i)).^2)./(2*vA));
g_Y_theta(1,i)=exp(sum(lg_Y_theta,1));%numerically more stable version
%g_Y_theta(1,i)=prod(exp(lg_Y_theta),1);

%% Step 3c
m(1,i)=(sum(lp_Y_theta,1)-sum(lg_Y_theta,1));
%m(1,i)=(log(p_Y_theta(1,i))-log(g_Y_theta(1,i)));%less stable version
end

avg_m=sum(m)/N;
for i=1:N
w(1,i)=avg_m * exp(m(1,i)-avg_m);
end
avg_w=sum(w)/N;
%% loglikelihood evaluation due to (Q)MLE
l_hat = llik_g + avg_m - log(N) + log(sum(exp(m(1,i)-avg_m)));
llik=mean(l_hat);
%l=  -(1/2)*T*log(2*pi) -(1/2)*sum(log(abs(mF)) +((vV.^2)./mF)); 
end
      