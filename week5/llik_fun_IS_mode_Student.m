function [llik]=llik_fun_IS_mode_Student(vy,theta,rand_eta,rand_u,N)
% performing the importance sampling for given parameters to obtain log-likelihood          
%% obtaining the mode - used for SPDK method as importance density
N=N;
sigma2Eta = theta(1);
phi = theta(3);
omega = theta(2);
v = theta(4);

T = length(vy);
I = ones(T,1);
g = 0;
q = zeros(T,1);
a0 = omega / (1 - phi);
P0 = sigma2Eta / (1 - phi^2);
c = 0;
d = omega;
mT = phi;
n=1;

%% MODE estimation
if g(n,1)==0 
     q = I + exp(-g).*(vy.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g - 0.5*vA + (q.^2);  
    Yn = z;
    H = A;        
    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2Eta,a0,P0);
    g = modsmooth;
    new_value=g(n);
     n=n+1;
else
end
while sqrt((new_value-g(n-1,1))^2)>0.00001 % for cycle is more suitable than while cycle; after 15 iterations convergence is reached surely  
        %while cycle
    q = I + exp(-g).*(vy.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g - 0.5*vA + (q.^2);  
    Yn = z;
    H = A;          
    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2Eta,a0,P0);
    g = modsmooth;
    new_value=g(n);
    n=n+1;
end

%% Step 1. Construct approximating g
  %based on chapter 14 in DK book about y_t close to 0
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
 [llik_g,~]=kf_smooth(y_star,A,mT,c,d,sigma2Eta,a0,P0);
%% Importance sampling and simulation smoother to obtain conditional mean of the signal and value of the log likelihood
vy_P = zeros(T,1);
vh_P = zeros(T+1,1);
vw = zeros(1,N);
vh_tilda = zeros(T,N);

% SPDK method for g(theta|Y_n), simulation smoother to simulate from g(theta|Y_n)  
for i = 1:N
% unconditional simulation - output: vh_P, vy_P 
h = omega /(1 - phi); %initialization value for vh since it is AR(1)
vh_P(1,1)= h; 
for t = 1:T
    eta = sqrt(sigma2Eta)*rand_eta(t,i);
    u = sqrt(vA(t))*rand_u(t,i);
    %vy_P(t,1) = exp(0.5*h)*u;
     vy_P(t,1) = h+c+u;
    h = omega + phi*h + eta;
    vh_P(t+1,1) = h;
end

% KFS to obtain vh_HP
[~,vh_HP] = kf_smooth(vy_P,A,phi,c,omega,sigma2Eta,a0,P0);

% KFS to obtain vtheta_H
vy_star = z;
[~,vh_H] = kf_smooth(vy_star,A,phi,c,omega,sigma2Eta,a0,P0);

% conditional simulation smoothing
vh_tilda(:,i) = vh_H - vh_HP + vh_P(1:end-1,1); %since the last value of vh_P is for t+1

%vh_tilda is a draw from importance density g(theta|Y_n)

%% Step 3b Evaluate p(Yn|theta) and g(Yn|theta)
%compute p(Yn|theta)
z_new=(vy-c)./(exp(0.5*omega/(1-phi)));
qt = 1 + exp(-vh_tilda(:,i)).*(z_new.^2)./((v-2));
lp_Y_theta=c-0.5.*(vh_tilda(:,i) + (v+1).*log(qt));
pdens =exp(sum(lp_Y_theta,1));%numerically more stable version
%p_Y_theta(1,i)=prod(exp(lp_Y_theta),1);

%compute g(Yn|theta)
lg_Y_theta=-0.5.*(log(2*pi*vA)) - (1./(2*vA)).*(Yn-c-vh_tilda(:,i)).^2;%not sure if c supposted to be here
%g_Ytheta=1/(sqrt(2*pi*vA)).*exp(-((y_star-theta_tilde(:,i)).^2)./(2*vA));
 gdens=exp(sum(lg_Y_theta,1));%numerically more stable version
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
mode = sum(vh_tilda.*exp(m-avg_m),2)./sum(exp(m-avg_m)); 
modeV = sum(((vh_tilda - mode*ones(1,N)).^2).* (ones(T,1)*exp(m-avg_m)),2)./sum(exp(m-avg_m)); 

%% loglikelihood evaluation due to (Q)MLE
l_hat = llik_g + avg_m - log(N) + log(sum(exp(m(1,i)-avg_m),2));
llik=mean(l_hat);
end