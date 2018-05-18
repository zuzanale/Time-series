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

print mean(y)

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
        lb = [10^(-7); 10^(-7); 10^(-7)]; 
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
plot1 = plot(h_smoothed);
set(plot1,'Color','r','LineWidth',1.5);
xticklabels({0,100,200,300,400,500,600,700,800,900})
axis([0 950 -20 5])


%% Ex. 4.1
%firstly for loop with 15 iterations (should be enough to converge to
%optimum
T = length(y);
I = ones(T,1);
g = zeros(T,1);

for i=1:15 %
    A = diag( (2 .*exp(g))./y.^2 );
    vA = (2 .*exp(g))./y .^2;
    z = g - (1/2).*vA+ I ;  
    Yn = z;
    H = A;    

[~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
g = modsmooth;

end  

figure(5)
subplot(2,1,1)
q2 = plot(h_smoothed);
set(q2,'Color','black','LineWidth',1.5)
hold on
q1 = plot(g);
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
q1 = plot(exp(0.5*g));
set(q1,'Color','black','LineWidth',2);
%axis([0 950 -20 5])
legend('Data','KF smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)

%% Ex. 4.1
% secondly the while loop, we used 0.00001 difference as a threshold
T = length(y);
I = ones(T,1);
g = zeros(T,1);
n=1;
if g(n,1)==0 
    A = diag( (2 .*exp(g))./y.^2 );
    vA = (2 .*exp(g))./y .^2;
    z = g - (1/2).*vA+ I ;  
    Yn = z;
    H = A;    
[~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
g = modsmooth;
else
    while sqrt(g(n,1)-g(n-1,1)^2)>0.00001 % for cycle is more suitable than while cycle; after 15 iterations convergence is reached surely  
        %while cycle
        A = diag( (2 .*exp(g))./y.^2 );
        vA = (2 .*exp(g))./y .^2;
        z = g - (1/2).*vA+ I ;  
        Yn = z;
        H = A;    

    [~,modsmooth] = kf_smooth(Yn,H,mT,c,d,sigma2_eta,a0,P0);
    g = modsmooth;
    n=n+1;
    end
end

figure(6)
subplot(2,1,1)
q2 = plot(h_smoothed);
set(q2,'Color','black','LineWidth',1.5)
hold on
q1 = plot(g);
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
q1 = plot(exp(0.5*g));
set(q1,'Color','black','LineWidth',2);
%axis([0 950 -20 5])
legend('Data','KF smoothed mean','Mode estimate');
xticklabels({0,100,200,300,400,500,600,700,800,900})
hold off
set(findall(gcf,'type','text'),'FontSize',14)
