%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 5: SV model (Assignment 3,4 included)
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

figure(1)
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
figure(3)
plot(x)
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold on
title('Visual representation of x(t)')
axis([0 950 -20 20])

mean(y)

%% Part C - QMLE
%% 2. Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
        %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
%% 3. Initial Parameter Values to parameters that has to be estimated
             
    sigma2_u = pi^2/2;  %this one is given to us, we do not have to estimate it this time Q = pi^2/2
    sigma2_eta_ini = 0.05;
    omega_ini = 0.05; 
    phi_ini = 0.8; 
    
    theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
        
%% 4. Parameter Space Bounds
        % parameter lower and upper bounds
        lb = [10^(-7); -100; 10^(-7)]; 
        ub = [inf; 10^7; 0.9999999];
      
      
%% 5. Optimize Log Likelihood Criterion

      [theta_hat,llik_val,exitflag]=...
          fmincon(@(theta) - llik_fun(x,theta),theta_ini,[],[],[],[],lb,ub,[],options);
    
display('parameter sigma2_eta')
theta_hat(1)
    
display('parameter omega')
theta_hat(2)
    
display('parameter phi')
theta_hat(3)

%% PART d
%% Our final estimates
  %  theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 3x1
    a0 = theta_hat(2)/(1-theta_hat(3));% for week 4 its zeta
    P0 = theta_hat(1)/(1- theta_hat(3)^2);
    c = -1.27; 
    H = (pi^2/2)*eye(size(data,1)); 
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

figure(4)
plot(x);
hold
plot(h_smoothed);
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 950 -20 5])

%% Ex. 4.1
% secondly the while loop, we used 0.00001 difference as a threshold
T = length(y);
I = ones(T,1);
g = zeros(T,1);
difference=10;
c=0;
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

figure(5)
subplot(2,1,1)
plot(h_smoothed);
hold on
plot(g);
legend('Kalman filter smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 900 -3 1])
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
plot(y);
hold on
plot(exp(0.5*h_smoothed));
hold on
plot(exp(0.5*g));
legend('Data','KF smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 900 -4 6])
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% Assign. 5.1 f) - the smoothed mean Ht based on the SV model for y_t
% by using importance sampling method.
%% 0. generate disturbances in order to have the same ones
    N = 100;
  rand_u = randn(T,N);  % generate a vector discrepances
  rand_eta = randn(T,N); %generate a vector discrepances 
%% QMLE
%% Estimating parameters sigma2u, omega, phi, sigma2eta via importance sampling

 %%  Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-9,... % argument convergence criteria
                         'MaxIter',10000); % maximum number of iterations    
        %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');

%% bounds for QML optimization
        % parameter lower and upper bounds
        lb = [10^(-7); -100; 10^(-7)]; 
        ub = [inf; 10^7; 0.9999999];

 %% first parameters from mode optimization from the last week 4
    
    sigma2_eta_ini = 0.0070;
    omega_ini =  -0.0070; 
    phi_ini =0.95;%0.9912;
    
    theta_ini = [sigma2_eta_ini;omega_ini; phi_ini]; %dimension 4x1
 %%mode calculation
[llik,mode,modeV]=llik_fun_IS(y,theta_ini,rand_u,rand_eta,N);

% Plot the results
figure(7)
subplot(2,1,1)
plot(g);
hold on
plot(mode);
hold on
plot(h_smoothed);
legend('Mode estimate','Mean estimate IS','Kalman filter smoothed mean');
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 900 -3 1])
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
plot(y);
hold on
plot(exp(0.5*h_smoothed));
hold on
plot(exp(0.5*g));
hold on
plot(exp(0.5*mode));
legend('Data','KF smoothed mean','Mode estimate','Mean with IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 900 -4 6])
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% Estimate the parameters by maximizing the simulated likelihood function
% for student disturbances
    sigma2_eta_ini = 0.0070;
    sigma_ini = 0.9; 
    phi_ini = 0.95;
    v_ini =5;
    N=100;
    theta_ini = [sigma2_eta_ini; sigma_ini; phi_ini; v_ini]; 

 %%  Optimization Options
      options = optimset('Display','iter',... %display iterations
                         'TolFun',1e-16,... % function value convergence criteria 
                         'TolX',1e-14,... % argument convergence criteria
                         'MaxIter',1000); % maximum number of iterations    
%% bounds for QML optimization
        % parameter lower and upper bounds
        lb = [10^(-7); 10^(-7); 10^(-7); 2.0001]; 
        ub = [10^7; 10^7; 0.999999999999; 300];
        
        
%% Optimize Log Likelihood Criterion

      [theta_hat,loglik,exitflag]=...
          fmincon(@(theta) -llik_fun_IS_Student_MLE(y,theta,rand_eta,rand_u,N),theta_ini,[],[],[],[],lb,ub,[],options);
      
[mode_st,modeV]=llik_fun_IS_Student(y,theta_ini,rand_eta,rand_u,N);
   
display('parameter sigma2_eta')
theta_hat(1)
    
display('parameter sigma')
theta_hat(2)
    
display('parameter phi')
theta_hat(3)

display('parameter nu')
theta_hat(4)

%plot the results
figure(8)
subplot(2,1,1)
plot(mode_st);
legend('Mean estimate IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

subplot(2,1,2)
plot(y);
hold on
plot(exp(0.5*mode_st));
legend('Data','Mean with IS');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)
