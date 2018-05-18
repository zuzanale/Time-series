function [llik]=llik_fun_IS_Student_MLE(vy,theta,rand_eta,rand_u,N)
% performing the importance sampling for given parameters to obtain log-likelihood          
%% obtaining the mode - used for SPDK method as importance density
sigma2_eta = theta(1);
phi = theta(3);
sigma= theta(2);
v = theta(4);

T = length(vy);
I = ones(T,1);
g = 0;
q = zeros(T,1);
a0 = 0;
P0 = sigma2_eta / (1 - phi^2);
c = 0;
d = 0;
mT = phi;
n=1;

%% MODE estimation
%% Step 1. Construct approximating g
  for t=1:T
      if abs(vy(t,1))<0.0000001
           vy(t,1)=0.001;
      else
          vy(t,1) = vy(t,1);
      end
  end
difference=10;  
vz=vy./sigma;
while difference>0.00001
%for i=1:15 % 1
     q = I + exp(-g).*(vz.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g - 0.5*vA + q;  
    Yn = z;
    H = A;        
    [llik_g,modsmooth] = kf_smooth_st(Yn,H,mT,0,0,sigma2_eta,a0,P0);
    difference=sqrt(sum((g-modsmooth).^2));
     g = modsmooth;
end
 y_star = z; 
%% Step 3 Construct log w
%% 3a) Draw theta^(i)
y_plus = zeros(T,1);
theta_plus= zeros(T+1,1);
theta_tilde = zeros(T,N);
p_Y_theta = zeros(1,N);
g_Y_theta = zeros(1,N);
m = zeros(1,N);

% SPDK method for g(theta|Y_n), simulation smoother to simulate from g(theta|Y_n)  
for i = 1:N

% KFS to obtain conditional theta_hat
[~,con_theta_hat] =kf_smooth_st(y_star,A,phi,0,0,sigma2_eta,0,P0);

% KFS to obtain uncoditional simulation theta^+
     alpha_plus(1,1)=a0;%initialize alpha^+
    theta_plus(1,1)=alpha_plus(1,1);%initialize theta+
    for t = 1:T
        u_plus(t,i) = sqrt(vA(t)).*rand_u(t,i);  % generate a vector discrepances
        eta_plus(t,i) = sqrt(sigma2_eta).*rand_eta(t,i); %generate a vector discrepances 
        alpha_plus(t+1,1)= phi * alpha_plus(t,1)+eta_plus(t,i);
        theta_plus(t+1,1)=alpha_plus(t+1,1);% i dont think this one is necessary
        y_plus(t,1)=  theta_plus(t,1)+c+u_plus(t,i);
    end
    
% the smoothed mean from the simulate linear gaussian model
[~,con_theta_hat_plus] =kf_smooth_st(y_plus,A,phi,0,0,sigma2_eta,0,P0);

% uconditional simulation smoothing
theta_tilde(:,i) = con_theta_hat + theta_plus(1:end-1,1) -con_theta_hat_plus; %since the last value of vh_P is for t+1

%% Step 3b Evaluate p(Yn|theta) and g(Yn|theta)
%compute p(Yn|theta)
vz=vy./sigma;
qt = 1 + exp(-theta_tilde(:,i)).*(vz.^2)./(v-2);
lp_Y_theta=log((gamma((v+1)/2))/gamma(v/2)) - 0.5*log((v-2)*sigma^2)-0.5.*(theta_tilde(:,i) + (v+1).*log(qt));
p_Y_theta(1,i)=exp(sum(lp_Y_theta,1));%numerically more stable version
%p_Y_theta(1,i)=prod(exp(lp_Y_theta),1);

%compute g(Yn|theta)
lg_Y_theta=-0.5.*(log(2*pi*vA)) - (1./(2*vA)).*(y_star-c-theta_tilde(:,i)).^2;
g_Y_theta(1,i)=exp(sum(lg_Y_theta,1));%numerically more stable version
%% Step 3c
m(1,i)=(sum(lp_Y_theta,1)-sum(lg_Y_theta,1));
%m(1,i)=(log(p_Y_theta(1,i))-log(g_Y_theta(1,i)));%less stable version
end
avg_m=sum(m)/N;
mode = sum(theta_tilde.*exp(m-avg_m),2)./sum(exp(m-avg_m)); 
modeV = sum(((theta_tilde - mode*ones(1,N)).^2).* (ones(T,1)*exp(m-avg_m)),2)./sum(exp(m-avg_m)); 

%% loglikelihood evaluation due to (Q)MLE
l_hat = llik_g + avg_m - log(N) + log(sum(exp(m-avg_m),2));
llik=mean(l_hat);
end