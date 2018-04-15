function [llik]=llik_fun_IS_mode(x,theta,rand_eta,rand_u,N)     
%debuggin
    %theta(1)=sigma2_eta_ini ;
    %theta(2)=omega_ini; 
    %theta(3)=phi_ini ;
    %x=y;
   % rand_eta=rand_eta;
   % rand_u=rand_u;
  
%theta_ini = [sigma2_eta_ini; omega_ini; phi_ini]; %dimension 3x1
a0=0;
        %a0=theta(2)/(1-theta(3));
        P0 = theta(1)/(1-theta(3)^2);
        %H =theta(1)*eye(size(x,1)); %this one might not be necessary to compute
        sigmaEta= theta(1);
        c = 0;%mean 
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
    %mH = H;
    mR = 1;

    %% initial values
    a = a0;
    p = P0;
    
%% MODE estimation
%standardized observations
difference=10;
while difference>0.00001
    %sqrt(sum(g-old_value).^2)>0.00001
    %vy_std=(vy-c)./(exp(0.5*omega/(1-phi)));
      A = diag( (2 .*exp(g))./vy.^2 );
      vA = (2 .*exp(g))./vy .^2;
      %z = g - (1/2).*vA+ I ;  
    z= g + I - exp(g)./vy.^2;
      Yn = z;
      H = A;    
    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
    difference=sqrt(sum((g-modsmooth).^2));
    g = modsmooth;
end
%% Step 1. Construct approximating g
  for t=1:T
      if z(t,1)<0.0000001
           z(t,1)=0.001;
      else
          z(t,1) = z(t,1);
      end
  end
  y_star = z;
  H = A;%this is not necessary, should be erased for consistency    
  
  %% Step 2. Evaluate log g(Yn)
theta_g=[sigma2_eta;omega;phi];
 %[llik_g=kf_smooth(y_star,A,mT,c,d,sigma2_eta,a0,P0); %evaluate log[g(Yn)] with H=A
llik_g=llik_fun_app(y_star,theta_g,H);     
%% Step 3 Construct log w
%% 3a) Draw theta^(i)
y_plus = zeros(T,1);
theta_plus= zeros(T+1,1);
theta_tilde = zeros(T,N);
p_Y_theta = zeros(1,N);
g_Y_theta = zeros(1,N);
m = zeros(1,N);
w = zeros(1,N);
a0=0;
P0 = theta(1)/(1-theta(3)^2);
for i = 1:N 

% KFS to obtain conditional theta_hat
[~,con_theta_hat] =kf_smooth(y_star,A,phi,c,omega,sigma2_eta,a0,P0);

% KFS to obtain uncoditional simulation theta^+
     %alpha_plus(1,1)=omega + phi *a0+sqrt(P0)*rand_eta(1,i);%initialize alpha^+
     alpha_plus(1,1)=a0;%initialize alpha^+
    theta_plus(1,1)=alpha_plus(1,1);%initialize theta+
    for t = 1:T
        u_plus(t,i) =  exp(0.5*(omega/(1-phi))).*exp(0.5*alpha_plus(t,1))*rand_u(t,i);  % generate a vector discrepances
        eta_plus(t,i) = sqrt(sigma2_eta).*rand_eta(t,i); %generate a vector discrepances 
        alpha_plus(t+1,1)= phi * alpha_plus(t,1)+eta_plus(t,i);
        theta_plus(t+1,1)=alpha_plus(t+1,1);% i dont think this one is necessary
        y_plus(t,1)=  theta_plus(t,1)+c+u_plus(t,i);
    end

% the smoothed mean from the simulate linear gaussian model
[~,con_theta_hat_plus] =kf_smooth(y_plus,A,phi,c,omega,sigma2_eta,a0,P0);

% uconditional simulation smoothing
theta_tilde(:,i) = con_theta_hat + theta_plus(1:end-1,1) -con_theta_hat_plus; %since the last value of vh_P is for t+1

%% Step 3b Evaluate p(Yn|theta) and g(Yn|theta)
%compute p(Yn|theta)
%z_new=(vy-c)./(exp(0.5*omega/(1-phi)));
lp_Y_theta=0.5.*(log(2*pi*exp(0.5*omega/(1-phi))^2) + theta_tilde(:,i) + vy.^2.*exp(-theta_tilde(:,i)));
p_Y_theta(1,i)=exp(sum(lp_Y_theta,1));%numerically more stable version
%p_Y_theta(1,i)=prod(exp(lp_Y_theta),1);

%compute g(Yn|theta)
%vA is always positive so it doesnt have to be in abs()
lg_Y_theta=-0.5.*(log(2*pi*vA)) - (1./(2*vA)).*(y_star-c-theta_tilde(:,i)).^2;
g_Y_theta(1,i)=exp(sum(lg_Y_theta,1));%numerically more stable version
%g_Y_theta(1,i)=prod(exp(lg_Y_theta),1);

%% Step 3c
m(1,i)=(sum(lp_Y_theta,1)-sum(lg_Y_theta,1));
%m(1,i)=(log(p_Y_theta(1,i))-log(g_Y_theta(1,i)));%less stable version
end

avg_m=sum(m)/N;
for i=1:N
w(1,i)=exp(avg_m) * exp(m(1,i)-avg_m);
end
avg_w=sum(w)/N;
mode = sum(theta_tilde.*exp(m-avg_m),2)./sum(exp(m-avg_m)); 
modeV = sum(((theta_tilde - mode*ones(1,N)).^2).* (ones(T,1)*exp(m-avg_m)),2)./sum(exp(m-avg_m)); 

%% loglikelihood evaluation due to (Q)MLE
l_hat = llik_g + avg_m - log(N) + log(sum(exp(m(1,i)-avg_m),2));
llik=mean(l_hat);
end
      