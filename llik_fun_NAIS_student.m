function [mode,modeV]=llik_fun_NAIS_student(vy,theta_hat,rand_eta,rand_u,N,M,threshold)     
%debugging
%theta_hat=theta;
%vy=y;

    [GH.z, GH.h]  = hernodes(M);
    T = size(vy,1);
    
    %initialize b and C
    b=zeros(T,1);
    mC=1*eye(size(vy,1)); % V = (Z*P1*Z')';
    C=ones(T,1);
    z =GH.z;      % the Gauss-Hermite abscissae z 
    h = GH.h;      % and the associated weights h(z)
    k=0;
    difference_b = 1;
    difference_C = 1;

%% Repeat until convergence

%% compute mean theta_hat and variance Vt
%standardized observations
while or((difference_b>=threshold),(difference_C>=threshold))
    old_b=b;
    old_C=C;
    k=k+1;
    y_star = b./C; 
    a0 = 0;
    P0 = theta_hat(1)/(1- theta_hat(3)^2);
    c = 0; 
    H =(mC).^(-1);
    sigma=theta_hat(2);
    phi=theta_hat(3);
    sigma2_eta=theta_hat(1);
    v=theta_hat(4);
    [h_smoothed,V] = kf_smooth_st_NAIS(y_star,H,phi,0,0,sigma2_eta,a0,P0);
    
    %% for each t
    %% compute theta^(t,j)
        theta_tilde = zeros(T,M);
        lp_Y_theta = zeros(1,M);
        for t=1:T
                theta_tilde(t,:) = h_smoothed(t) + V(t).^(1/2).*z(:);
                vz(t)=vy(t)/sigma;
                qt = 1 + exp(-theta_tilde(t,:)).*(vz(t).^2)./(v-2);
                lp_Y_theta(t,:)=log((gamma((v+1)/2))/gamma(v/2)) - 0.5*log((v-2)*sigma^2)-0.5.*(theta_tilde(t,:) + (v+1).*log(qt));
                %lp_Y_theta(t,:)=log(gamma((v+1)/2)) - 0.5*log(v-2) - 0.5*log(pi) - log(gamma(v/2))-0.5.*(theta_tilde(t,:) + (v+1).*log(qt));
                %lp_Y_theta(t,:)=-0.5.*(log(2*pi*exp(0.5*omega/(1-phi))^2) + theta_tilde(t,:) + vy(t).^2.*exp(-theta_tilde(t,:)));
            %weigheted least squares regression
            x=[1*ones(1,M);theta_tilde(t,:);-0.5.*theta_tilde(t,:).^2]';
            W = diag(exp(0.5.*(z.^2)).* h)';
            beta(t,:)=inv(x'*W*x)*x'*W*lp_Y_theta(t,:)';
            b(t)=beta(t,2);
            C(t)=beta(t,3);
        end
        mC=diag(C);
        difference_b=mean((old_b-b).^2);
        difference_C=mean((old_C-C).^2);
 end
 
%% Step 3c
%% 3a) Draw theta^(i)
y_star = b./C; 
y_plus = zeros(T,1);
theta_plus= zeros(T+1,1);
theta_tilde2 = zeros(T,N);
p_Y_theta = zeros(1,N);
g_Y_theta = zeros(1,N);
log_w = zeros(1,N);
w= zeros(1,N);
for i = 1:N 

% KFS to obtain conditional theta_hat
[~,con_theta_hat] = kf_smooth_st(y_star,(mC)^(-1),phi,0,0,sigma2_eta,0,P0);

% KFS to obtain uncoditional simulation theta^+
    alpha_plus(1,1)=a0;%initialize alpha^+
    theta_plus(1,1)=alpha_plus(1,1);%initialize theta+
    for t = 1:T
        u_plus(t,i) = sqrt((C(t))^(-1)).*rand_u(t,i);  % generate a vector discrepances
        eta_plus(t,i) = sqrt(sigma2_eta).*rand_eta(t,i); %generate a vector discrepances 
        alpha_plus(t+1,1)= phi * alpha_plus(t,1)+eta_plus(t,i);
        theta_plus(t+1,1)=alpha_plus(t+1,1);% i dont think this one is necessary
        y_plus(t,1)=  theta_plus(t,1)+u_plus(t,i);
    end
% the smoothed mean from the simulate linear gaussian model
[~,con_theta_hat_plus] =kf_smooth_st(y_plus,(mC)^(-1),phi,0,0,sigma2_eta,0,P0);

% uconditional simulation smoothing
theta_tilde2(:,i) = con_theta_hat + theta_plus(1:end-1,1) -con_theta_hat_plus; %since the last value of vh_P is for t+1

%% Step 3b Evaluate p(Yn|theta) and g(Yn|theta)
%compute p(Yn|theta)
vz=vy./sigma;
qt = 1 + exp(-theta_tilde2(:,i)).*(vz.^2)./(v-2);
lp_Y_theta=log((gamma((v+1)/2))/gamma(v/2)) - 0.5*log((v-2)*sigma^2)-0.5.*(theta_tilde2(:,i) + (v+1).*log(qt));
%lp_Y_theta=log(gamma((v+1)/2)) - 0.5*log(v-2) - 0.5*log(pi) - log(gamma(v/2))-0.5.*(theta_tilde(:,i) + (v+1).*log(qt));     
p_Y_theta(1,i)=exp(sum(lp_Y_theta,1));%numerically more stable version

%compute g(Yn|theta)
a=0.5.*(log(abs(C))-log(2*pi)-b.*y_star);
lg_Y_theta=a + b.*theta_tilde2(:,i) -0.5.*(theta_tilde2(:,i).^2).*C;
g_Y_theta(1,i)=exp(sum(lg_Y_theta,1));%numerically more stable version

%% Step 3c
log_w(1,i)=(sum(lp_Y_theta,1)-sum(lg_Y_theta,1));
w(1,i)=p_Y_theta(1,i)./g_Y_theta(1,i);
%m(1,i)=(log(p_Y_theta(1,i))-log(g_Y_theta(1,i)));%less stable version
end
m=log_w;
avg_m=sum(m)/N;
mode = sum(theta_tilde2.*exp(m-avg_m),2)./sum(exp(m-avg_m)); 
modeV = sum(((theta_tilde2 - mode*ones(1,N)).^2).* (ones(T,1)*exp(m-avg_m)),2)./sum(exp(m-avg_m)); 
%for j=1:M
%    theta_tilde(:,j)
%    g_theta_y=exp(-0.5.*(z(j).^2))./(sqrt(2*pi.*V));
%    m(1,j)=h(j).*exp(z(j).^2);
%m(1,j)=(sum(lp_Y_theta,1)-sum(lg_Y_theta,1));
%m(1,i)=(log(p_Y_theta(1,i))-log(g_Y_theta(1,i)));%less stable version
%end
end
      