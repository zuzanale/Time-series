%%  TIME SERIES ECONOMETRICS
%
%   Final assignment: SV model
%
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1)
%   April 2018

%% 0. Clean Workspace and Command Window

clear all        %clear workspace
clc              %clear command window

%% 0. Read Data
data = xlsread('SPData.xlsx');

%% 1a: Transformation of level into returns

r = data(:,3);
log_r = log(r);
diff_r = diff(log_r);
log_y = 100*diff_r;
y = log_y - mean(log_y);
%compute the realized measures
RM = data(:,2);

% Convert Excel to matlab dates
Dates = data(:,1) + 693960;
dates = datestr(Dates);


% plots of the prices, log prices and log returns
figure(1)
subplot(3,1,1)
plot(Dates,r)
datetick('x',10)
title('Prices')
hold on
subplot(3,1,2)
plot(Dates,log_r)
datetick('x',10) 
title('Log Prices')
hold on
subplot(3,1,3)
plot(Dates(2:end,1),y)
datetick('x',10)
title('Log Returns')
hold off

% descriptive statistics of the prices, log prices and log returns
series = {'Price';'Log price';'Log returns'};
Mean = [mean(r); mean(log_r); mean(y)];
Median = [median(r); median(log_r); median(y)];
Min = [min(r); min(log_r); min(y)];
Max = [max(r); max(log_r); max(y)];
SD = [std(r); std(log_r); std(y)];
Skewness = [skewness(r); skewness(log_r); skewness(y)];
ExKurtosis = [kurtosis(r)-3; kurtosis(log_r)-3; kurtosis(y)-3];

f = figure('Position',[100 100 400 250]);
dat = [{'Price';'Log price';'Log returns'} num2cell(Mean) num2cell(Median) num2cell(Min) num2cell(Max) num2cell(SD) num2cell(Skewness) num2cell(ExKurtosis)];
colname =   {' ', 'Mean', 'Median', 'Min', 'Max', 'SD', 'Skewness', 'ExKurtosis'};
colformat = {'char', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
t = uitable('Units','normalized','Position',[0.1 0.1 0.7 0.8],'Data', dat,... 
  'RowName',{'1','2','3'}, 'ColumnName',colname,'ColumnFormat', colformat);

%% 1b: Linear SV model
x = log(y.^2); %demeaning incorporated

% plot x
figure(4)
plot(Dates(2:end,1),x)
datetick('x',10)
hold on
title('Visual representation of x(t)')

%% 1c: Estimating coefficients by QML methods
%% 1. Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
        %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
%% 2. Initial Parameter Values to parameters that has to be estimated
             
    sigma2_u = pi^2/2;  %this one is given to us, we do not have to estimate it this time
    sigma2_eta_ini = 0.05;
    omega_ini = 0.05; 
    phi_ini = 0.8; 
    
    theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
        
%% 3. Parameter Space Bounds
        % parameter lower and upper bounds
        lb = [10^(-7); -100; 10^(-7)]; 
        ub = [inf; 10^7; 0.9999999];
      
      
%% 4. Optimize Log Likelihood Criterion

      [theta_hat,llik_val,exitflag]=...
          fmincon(@(theta) - llik_fun(x,theta),theta_ini,[],[],[],[],lb,ub,[],options);
    
display('parameter sigma2_eta')
theta_hat(1)
    
display('parameter omega')
theta_hat(2)
    
display('parameter phi')
theta_hat(3)

%% 1d: Computing the smoothed mean
%% Our final estimates
  %  theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
    a0 = theta_hat(2)/(1-theta_hat(3));
    P0 = theta_hat(1)/(1- theta_hat(3)^2);
    c = -1.27; 
    H =(pi^2/2)*eye(size(y,1));
    d = theta_hat(2);
    mT = theta_hat(3);
    omega=theta_hat(2);
    phi=theta_hat(3);
    sigma2_eta=theta_hat(1);
 
    [llik,h_smoothed] = kf_smooth(x,H,phi,c,omega,sigma2_eta,a0,P0);

    vresid=NaN(length(x),1);
    vresid = x - h_smoothed;
    vresid_mean=mean(vresid);%this is the same as the mean of disturbances 
                               %u_t 

est = {'sigma2eta';'omega';'phi'};
Estimate=[sigma2_eta; omega; phi];
tableEst = table(Estimate,'RowNames',est);

figure(5)
plot(Dates(2:end,1),x)
datetick('x',10)
hold on
plot(Dates(2:end,1),h_smoothed);
hold off


%% 1e: Computing the smoothed mode
% the while loop, we used 0.00001 difference as a threshold
T = length(y);
I = ones(T,1);
g = zeros(T,1);
difference=10;
a0=0;
while difference>0.00001*ones(T,1) %for each time t the difference should be smaller than our threshold
    %y_std=(vy-c)./(exp(0.5*omega/(1-phi)));
      A = diag( (2 .*exp(g))./y.^2 );
      vA = (2 .*exp(g))./y .^2;
        z= g + I - exp(g)./y.^2;
      Yn = z;
      H = A;    
    [~,modsmooth] = kf_smooth(Yn,H,mT,0,0,sigma2_eta,0,P0);
    difference=sqrt(((g-modsmooth).^2));
    g = modsmooth;
end

figure(6)
subplot(2,1,1)
plot(Dates(2:end,1),h_smoothed);
datetick('x',10)
hold on
plot(Dates(2:end,1),g);
legend('Kalman filter smoothed mean','Mode estimate');
hold off
set(findall(gcf,'type','text'),'FontSize',14)


subplot(2,1,2)
plot(Dates(2:end,1),y);
datetick('x',10)
hold on
plot(Dates(2:end,1),exp(0.5*h_smoothed));
hold on
plot(Dates(2:end,1),exp(0.5*g));
legend('Data','KF smoothed mean','Mode estimate');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% 1f. Smoothed mean by using SPDK and NAIS method
%% SPDK method 
%should also simulated likelihood function there????
 % by using importance sampling method.
%% 0. generate disturbances in order to have the same ones
    N = 100;
  rand_u = randn(T,N);  % generate a vector discrepances
  rand_eta = randn(T,N); %generate a vector discrepances 

 %% first parameters from mode optimization from the last week 4
    
    %sigma2_eta_ini = 0.0070;
    %omega_ini =  -0.0070; 
    %phi_ini =0.9912;
    theta_ini = theta_hat; % use our QML estimates from c
[llik,mode,modeV]=llik_fun_IS(y,theta_ini,rand_u,rand_eta,N);

% compute mode with NAIS
M=30;
N=100;
threshold = 0.0001; % convergence tolerance
[mode_nais,modeV_nais]=llik_fun_NAIS(y,theta_ini,rand_eta,rand_u,N,M,threshold);          

 %% 1g.Compare the three different estimates for h_t(display them in graphs)    
figure(7)
subplot(2,1,1)
plot(Dates(2:end,1),g);
datetick('x',10)
hold on
plot(Dates(2:end,1),mode);
hold on
plot(Dates(2:end,1),h_smoothed);
hold on
plot(Dates(2:end,1),mode_nais);
legend('Mode estimate','Mean estimate IS','Kalman filter smoothed mean','Mean estimate NAIS');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
plot(Dates(2:end,1),y);
datetick('x',10)
hold on
plot(Dates(2:end,1),exp(0.5*h_smoothed));
hold on
plot(Dates(2:end,1),exp(0.5*g));
hold on
plot(Dates(2:end,1),exp(0.5*mode));
hold on
plot(Dates(2:end,1),exp(0.5*mode_nais));
legend('Data','KF smoothed mean','Mode estimate','Mean with IS','Mean with NAIS');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% 1.h How would you modify the approximate maximum likelihood method for estimating the unknown
%parameters?
%% 1.i Estimate the parameters by maximizing the simulated likelihood function for the Student's t SV
%model using an importance sampling method.

%% Estimating parameters sigma2u, omega, phi, v, sigma2eta via importance sampling
% initial values          
   sigma2_eta_ini = 0.0070;
    sigma_ini = 0.9; 
    phi_ini = 0.95;
    v_ini =5;
    N=100;
    theta_ini = [sigma2_eta_ini; sigma_ini; phi_ini; v_ini]; 

 %%  Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
%% bounds for QML optimization
        % parameter lower and upper bounds
       lb = [10^(-7); 10^(-7); 10^(-7); 2.0001]; 
       ub = [10^7; 10^7; 0.999999999999; 300];
        
%% Optimize Log Likelihood Criterion

      [theta_hat,loglik,exitflag]=...
          fmincon(@(theta) -llik_fun_IS_Student_MLE(y,theta,rand_eta,rand_u,N),theta_ini,[],[],[],[],lb,ub,[],options);
    
display('parameter sigma2_eta')
sigma2_eta=theta_hat(1)
    
display('parameter omega')
sigma=theta_hat(2)
    
display('parameter phi')
phi=theta_hat(3)

display('parameter nu')
v=theta_hat(4)

%% 1.j Redo questions (e) - (f) but now for the Student's t
%SV model and the ML estimates of (i).

%(e) smoothed mode
T = length(y);
I = ones(T,1);
g_st = zeros(T,1);
difference=10;
a0=0;
vz=y/sigma;%standartized observations
while difference>0.00001
     q = I + exp(-g_st).*(vz.^2)./((v-2).*I);
    A = diag(2.*I ./((v+1).*(q - I)).*(q.^2));
    vA =  2.*I./((v+1).*(q - I)).*(q.^2);
    z = g_st - 0.5*vA + q;  
    Yn = z;
    H = A;        
    [~,modsmooth] = kf_smooth_st(Yn,H,mT,0,0,sigma2_eta,a0,P0);
    difference=sqrt(sum((g_st-modsmooth).^2));
     g_st = modsmooth;
end

% 1f. Smoothed mean by using SPDK and NAIS method
% SPDK method 
[mode_st,modeV_st]=llik_fun_IS_Student(y,theta_hat,rand_eta,rand_u,N);

% NAIS method 
M=30;
N=100;
threshold = 0.01; % convergence tolerance
[mode_nais_st,modeV_nais_st]=llik_fun_NAIS_student(y,theta_hat,rand_eta,rand_u,N,M,threshold);          

figure(8)
subplot(2,1,1)
title('Student density model')
plot(Dates(2:end,1),modsmooth);
datetick('x',10)
hold on
plot(Dates(2:end,1),mode_st);
hold on
plot(Dates(2:end,1),mode_nais_st);
legend('Mode estimate','Mean estimate IS','Mean estimate NAIS');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
p1=plot(Dates(2:end,1),y);
datetick('x',10)
hold on
p2=plot(Dates(2:end,1),exp(0.5*modsmooth));
hold on
p3=plot(Dates(2:end,1),exp(0.5*mode_st));
hold on
p4=plot(Dates(2:end,1),exp(0.5*mode_nais_st));
legend('Data','Mode estimate','Mean with IS','Mean with NAIS');
hold off
set(p1,'LineWidth',0.1);
set(p2,'LineWidth',1.3);
set(p3,'LineWidth',1.3);
set(p4,'LineWidth',1.3);
set(findall(gcf,'type','text'),'FontSize',14)

%% 1.k Compare your estimates and analysis between the Gaussian and Student's t SV models.
figure(9)
title('Normal density model')
subplot(2,2,1)
plot(Dates(2:end,1),g);
datetick('x',10)
hold on
plot(Dates(2:end,1),mode);
hold on
plot(Dates(2:end,1),h_smoothed);
hold on
plot(Dates(2:end,1),mode_nais);
legend('Mode estimate','Mean estimate IS','Kalman filter smoothed mean','Mean estimate NAIS');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,2,2)
q1=plot(Dates(2:end,1),y);
datetick('x',10)
hold on
q2=plot(Dates(2:end,1),exp(0.5*h_smoothed));
hold on
q3=plot(Dates(2:end,1),exp(0.5*g));
hold on
q4=plot(Dates(2:end,1),exp(0.5*mode));
hold on
q5=plot(Dates(2:end,1),exp(0.5*mode_nais));
legend('Data','KF smoothed mean','Mode estimate','Mean with IS','Mean with NAIS');
hold off
set(q1,'LineWidth',0.1);
set(q2,'LineWidth',1.3);
set(q3,'LineWidth',1.3);
set(q4,'LineWidth',1.3);
set(q5,'LineWidth',1.3);
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,2,3)
title('Student density model')
plot(Dates(2:end,1),modsmooth);
datetick('x',10)
hold on
plot(Dates(2:end,1),mode_st);
hold on
plot(Dates(2:end,1),mode_nais_st);
legend('Mode estimate','Mean estimate IS','Mean estimate NAIS');
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,2,4)
p1=plot(Dates(2:end,1),y);
datetick('x',10)
hold on
p2=plot(Dates(2:end,1),exp(0.5*modsmooth));
hold on
p3=plot(Dates(2:end,1),exp(0.5*mode_st));
hold on
p4=plot(Dates(2:end,1),exp(0.5*mode_nais_st));
legend('Data','Mode estimate','Mean with IS','Mean with NAIS');
hold off
set(p1,'LineWidth',0.1);
set(p2,'LineWidth',1.3);
set(p3,'LineWidth',1.3);
set(p4,'LineWidth',1.3);
set(findall(gcf,'type','text'),'FontSize',14)
