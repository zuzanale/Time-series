
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

%% returns transformation
vy = data; %demeaned data

Dates = data(: ,1) + datenum ('30 DEC1899'); 
dates = datestr(Dates);

n = 1;
T=size(data,1);
N=100;
par = zeros(5,n);
RN_u = randn(T,N);  % generate a vector discrepances
RN_eta = randn(T,N); %generate a vector discrepances 

%% Estimating parameters sigma2u, omega, phi, v, sigma2eta via importance sampling
% initial values          
    sigma2_u_ini = 2;
    sigma2_eta_ini = 0.001;
    omega_ini = 0.0002; 
    phi_ini = 0.99;
    v_ini = 30;
    N=100;
    theta_ini = [sigma2_eta_ini; omega_ini; phi_ini; v_ini]; 

 %%  Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
%% bounds for QML optimization
        % parameter lower and upper bounds
        lb = [10^(-10); 10^(-10); 10^(-10); 2.0001]; 
        ub = [10^10; 10^9; 0.9999999; 100];
        
%% Optimize Log Likelihood Criterion

      [vparameter_est,loglik,exitflag]=...
          fmincon(@(vparameter) -llik_fun_IS_mode_Student(vy,vparameter,RN_eta,RN_u,N),theta_ini,[],[],[],[],lb,ub,[],options);

% transformation of final parameters
sigma2_eta_est=vparameter_est(1);
omega_est = vparameter_est(2);
phi_est = vparameter_est(3); 
v_est =  vparameter_est(4);

%mode estimation
[llik,mode,modeV]=llik_fun_IS_Student(vy,vparameter_est,RN_eta,RN_u,N);
% table with estimates
est = {'sigma2eta';'omega';'phi';'v'};
Estimate=[sigma2_eta_est; omega_est; phi_est; v_est];
estTable = table(Estimate,'RowNames',est);

figure(8)
subplot(2,1,1)
q1 = plot(mode);
set(q1,'Color','blue','LineWidth',1.5)
legend('Mode estimate IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
qd = plot(vy);
set(qd,'Color','black')
hold on
%q2 = plot(exp(0.5*h_smoothed));
%set(q2,'Color','red','LineWidth',2)
%hold on
%q1 = plot(exp(0.5*g));
%set(q1,'Color','blue','LineWidth',2);
%hold on
q3 = plot(exp(0.5*mode));
set(q3,'Color','blue','LineWidth',2);
%axis([0 950 -20 5])
legend('Data','Mode with IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% PARTS A)-E)
%% Part a) - plot of the given returns
y = data;

figure(9)
y1 = subplot(2,1,1)
plot(y)
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold on
title('Visual representation of the returns')
y2 = subplot(2,1,2)
histogram(y)
hold off
axis([y1],[0 950 -5 5])
axis([y2],[-5 5 0 200])
title('Histogram of the returns')

series = {'Returns'};
Mean = [ mean(y)];
Median = [ median(y)];
Min = [ min(y)];
Max = [ max(y)];
SD = [std(y)];
Skewness = [ skewness(y)];
ExKurtosis = [ kurtosis(y)-3];

f = figure('Position',[100 100 400 250]);
dat = [{'Returns'} num2cell(Mean) num2cell(Median) num2cell(Min) num2cell(Max) num2cell(SD) num2cell(Skewness) num2cell(ExKurtosis)];
colname =   {' ', 'Mean', 'Median', 'Min', 'Max', 'SD', 'Skewness', 'ExKurtosis'};
colformat = {'char', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
t = uitable('Units','normalized','Position',[0.1 0.1 0.7 0.8],'Data', dat,... 
  'RowName',{'1'}, 'ColumnName',colname,'ColumnFormat', colformat);

%% Part b)
x = log((y-mean(y)).^2); %demeaning incorporated

% plot x
figure(10)
plot(x)
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold on
title('Visual representation of x(t)')
axis([0 950 -20 20])

mean(y)

%% PART d
%% Our final estimates
  %  theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
    a0 = theta_hat(2)/(1-theta_hat(3));% for week 4 its zeta
    P0 = theta_hat(1)/(1- theta_hat(3)^2);
    c = 0; 
    H = (pi^2/2)*eye(size(data,1)); 
    d = theta_hat(2); 
    mT = theta_hat(3);
    omega=omega_est ;
    phi=phi_est;
    sigma2_eta=sigma2_eta_est;
 
    [llik,h_smoothed] = kf_smooth(x,H,phi,c,omega,sigma2_eta,a0,P0);

    vresid=NaN(length(x),1);
    vresid = x - h_smoothed;
    vresid_mean=mean(vresid);%this is the same as the mean of disturbances 
                               %u_t 

est = {'sigma2eta';'omega';'phi'};
Estimate=[sigma2_eta; omega; phi];
tableEst = table(Estimate,'RowNames',est);

figure(11)
plot(x);
hold
plot1 = plot(h_smoothed);
set(plot1,'Color','r','LineWidth',1.5);
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 950 -20 5])


%% Ex. 4.1
%firstly for loop with 15 iterations (should be enough to converge to
%optimum
T = length(y);
I = ones(T,1);
g = 0;
q = zeros(T,1);
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

figure(12)
subplot(2,1,1)
q2 = plot(g);
set(q2,'Color','black','LineWidth',1.5)
hold on
q1 = plot(mode);
set(q1,'Color','blue','LineWidth',1.5)
hold on
q3 = plot(h_smoothed);
set(q3,'Color','red','LineWidth',1.5)
legend('Mode estimate','Mode estimate IS','Kalman filter smoothed mean');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
qd = plot(y);
set(qd,'Color','grey')
hold on
q2 = plot(exp(0.5*h_smoothed));
set(q2,'Color','red','LineWidth',2)
hold on
q1 = plot(exp(0.5*g));
set(q1,'Color','blue','LineWidth',2);
hold on
q3 = plot(exp(0.5*mode));
set(q3,'Color','black','LineWidth',2);
%axis([0 950 -20 5])
legend('Data','KF smoothed mean','Mode estimate','Mode with IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)
