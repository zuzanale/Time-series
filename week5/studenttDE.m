%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 4: SV model (Assignment 3 included)
%
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1)
%   March 2018 

%% 0. Clean Workspace and Command Window

clear all        %clear workspace
clc              %clear command window

%% 0. Read Data
A=importdata('sv.dat');
A.data;
fid = fopen('sv.dat','r');
datacell = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 1, 'Collect', 1);
fclose(fid);
A.data = datacell{1};
data=A.data(:,1);

%% Part a) - plot of the given returns
y = data;


%% Part b)
x = log((y-mean(y)).^2); %demeaning incorporated



%% Part C - QMLE
%% 2. Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
        %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
%% 3. Initial Parameter Values to parameters that has to be estimated
             
    theta_hat(1) = 0.0011;
    theta_hat(2)= 0.0032; 
    theta_hat(3) = 0.9917; 
    theta_hat(4) = 36.0892;
    
%% PART d
%% Our final estimates
  %  theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
    a0 = theta_hat(2)/(1-theta_hat(3));% for week 4 its zeta
    P0 = theta_hat(1)/(1- theta_hat(3)^2);
    c = 0; 
    H = eye(size(data,1)); 
    d = theta_hat(2); 
    mT = theta_hat(3);
    omega=theta_hat(2);
    phi=theta_hat(3);
    sigma2_eta=theta_hat(1);
    nu = theta_hat(4);
 
    [llik,h_smoothed] = kf_smooth_studentt(x,H,phi,c,omega,sigma2_eta,a0,P0,nu);

    vresid=NaN(length(x),1);
    vresid = x - h_smoothed;
    vresid_mean=mean(vresid);%this is the same as the mean of disturbances 
                               %u_t 

est = {'sigma2eta';'omega';'phi'};
Estimate=[sigma2_eta; omega; phi];
tableEst = table(Estimate,'RowNames',est);

figure(4)
plot(x);
hold
plot1 = plot(h_smoothed);
set(plot1,'Color','r','LineWidth',1.5);
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 950 -20 5])


%% Ex. 4.1
% secondly the while loop, we used 0.00001 difference as a threshold
T = length(y);
I = ones(T,1);
g = 0;
q = zeros(T,1);
n=1;
v=theta_hat(4);
%% MODE estimation
if g(n,1)==0 
     q = I + exp(-g).*(y.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g - 0.5*vA + (q.^2);  
    Yn = z;
    H = A;        
    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
    g = modsmooth;
    new_value=g(n);
     n=n+1;
else
end
while sqrt((new_value-g(n-1,1))^2)>0.00001 % for cycle is more suitable than while cycle; after 15 iterations convergence is reached surely  
        %while cycle
    q = I + exp(-g).*(y.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g - 0.5*vA + (q.^2);  
    Yn = z;
    H = A;          
    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
    g = modsmooth;
    new_value=g(n);
    n=n+1;
end

figure(6)
subplot(2,1,1)
q2 = plot(h_smoothed);
set(q2,'Color','black','LineWidth',1.5)
hold on
q1 = plot(modsmooth);
set(q1,'Color','blue','LineWidth',1.5)
legend('Kalman filter smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
qd = plot(y);
set(qd,'Color','blue')
hold on
q2 = plot(exp(0.5*h_smoothed));
set(q2,'Color','red','LineWidth',2)
hold on
q1 = plot(exp(0.5*modsmooth));
set(q1,'Color','black','LineWidth',2);
%axis([0 950 -20 5])
legend('Data','KF smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)
