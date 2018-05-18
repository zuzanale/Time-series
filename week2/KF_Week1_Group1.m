%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 1: LOCAL LEVEL MODEL AND KALMAN FILTER
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1) 

%% 0. Clean Workspace and Command Window

clear all        %clear workspace
clc              %clear command window

%% 0. Read Data
A=importdata('Nile.dat');
A.data;
fid = fopen('Nile.dat','r');
datacell = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'HeaderLines', 1, 'Collect', 1);
fclose(fid);
A.data = datacell{1};
data=A.data(:,1);

%% 1. Setup

    T=size(data,1);  % sample size
    
%% PART A: Implementing Kalman filter for LLM
%% 2. Parameter Initialisation
    
    q = 0.1; %1469.1 for (d) and keep same for (f) (0.1)
    sigma_eps = 1; %15099 for (d) (1 for initial a)
    sigma_eta = q;
    p = 10^7;
    y = data;
    a = 0; %use y(1) for initialising at t=2 

%% 3. Generate Innovations
        
    epsilon = sqrt(sigma_eps)*randn(T,1);  % generate a vector of T random normal 
                                     % variables with variance sigma_eps^2  
       
    eta = sqrt(sigma_eta)*randn(T,1);      % generate a vector of T random normal variables
                                     % with variance sigma_eta^2                                
%% 4. Apply Kalman Filter
    
for t=1:T
   mu(t) = a;  
   p1(t) = p;
   v = y(t) - a;
   k = p / (p + sigma_eps);
   a = a + k*v;
   p = k*sigma_eps + sigma_eta;
   vt(t) = v;
   kt(t) = k;
end

%% 5. Plot Slide 65
x =1871:1:1970;
figure(1)
plot(x(2:100),mu(2:100),'b')  % plot the time-series y in blue 'b'  
hold on
plot(x(2:100),y(2:100),'r-o')

%% PART B: Replicate figure on slide 67
%% observation weights local level
 weight=zeros(21,1);
 prod_k=cumprod(1-kt(80:100));
 weight(21)=kt(79)*(1-kt(80));
for i=1:21
    weight(22-i)= kt(79)*prod_k(i);
end

%% Observation weights global level 
%global level model with fixed mu, weigthed sum of past observations   
%% Kalman filter    
for t=1:T
   mu=a;
   p1(t) = p;
   v = y(t) - a;
   k = p / (p + sigma_eps);
   a = a + k*v;
   p = k*sigma_eps;
   vt(t) = v;
   kt(t) = k;
end

%% Observation weights local level
globalweight=zeros(21,1);
 prod_k=cumprod(1-kt(80:100));
 globalweight(21)=kt(79)*(1-kt(80));
for i=1:21
    globalweight(22-i)= kt(79)*prod_k(i);
end

% plot graph with global and local weights
x =-20:1:0; % set limits for axis x
figure(2)
subplot(2,1,2)
bar(x,globalweight,'b') %plot global level
axis([-20 0 0 0.05])
title('Global Level')

subplot(2,1,1)
bar(x,weight,'b') %plot local level
axis([-20 0 0 0.3])
ylim([0.00 0.3])
title('Local Level: Filtered Weights')

%% PART C: different weight functions

    q = [10,1,0.1,0.001]; %1469.1 for (d) and keep same for (f) (0.1)
    weights=zeros(21,size(q,2));
    sigma_eps = 1; %15099 for (d) (1 for initial a)
    
%% Use kalman filter for each q

    for l=1:size(q,2)
        sigma_eta = q(l);
        p = 10^7;
        y = data;
        a = 0; %use y(1) for initialising at t=2 
        
        %% Apply Kalman Filter

        for t=1:T
           mu(t) = a;  
           p1(t) = p;
           v = y(t) - a;
           k = p / (p + sigma_eps);
           a = a + k*v;
           p = k*sigma_eps + sigma_eta;
           vt(t) = v;
           kt(t) = k;
        end
     prod_k=cumprod(1-kt(80:100));
     weights(21,l)=kt(79)*(1-kt(80));
        for i=1:21
            weights(22-i,l)= kt(79)*prod_k(i);
        end
    end
    x =-20:1:0; %set the axis x
    
    %plot data in 4 figures
    figure(3)
    
    subplot(2,2,1)
    bar(x,weights(:,1),'b')
    axis([-20 0 0 0.3])
    ylim([0.00 0.3])
    title('Local Level: Filtered Weights')
    xlabel('q = 10')
    
    subplot(2,2,2)
    bar(x,weights(:,2),'b')
    axis([-20 0 0 0.3])
    ylim([0.00 0.3])
    title('Local Level: Filtered Weights')
    xlabel('q = 1')
      
    subplot(2,2,3)
    bar(x,weights(:,3),'b')
    axis([-20 0 0 0.3])
    ylim([0.00 0.3])
    title('Local Level: Filtered Weights')
    xlabel('q = 0.1')
    
    subplot(2,2,4)
    bar(x,weights(:,4),'b')
    axis([-20 0 0 0.3])
    ylim([0.00 0.3])
    title('Local Level: Filtered Weights')
    xlabel('q = 0.001')
    
%% PART D: Apply Kalman filter to Nile Data

%% Parameter Initialisation
    
    q = 1469.1 ;
    sigma_eps = 15099; 
    sigma_eta = q;
    p = 10^7;
    a = 0; 
                      
%% Apply Kalman Filter
    
for t=1:T
   mu(t) = a;  
   p1(t) = p;
   v = y(t) - a;
   k = p / (p + sigma_eps);
   a = a + k*v;
   p = k*sigma_eps + sigma_eta;
   vt(t) = v;
   kt(t) = k;
end

%% PART E: Replicate Slide 65
x =1871:1:1970;
figure(4)
plot(x(2:100),mu(2:100),'b')  % plot the time-series y in blue 'b'  
hold on
plot(x(2:100),y(2:100),'r-o')

%% PART F: Initialize at t=2
    %q = 1469.1 ;   %stays the same as in the previous part
    %sigma_eps = 15099;  %stays the same as in the previous part
    %sigma_eta = q;     %stays the same as in the previous part
    p2 = sigma_eps+sigma_eta;
    a2 = y(1); 
                            
%% Apply Kalman Filter
    
for t=2:T
   mu2(t) = a2;  
   p1(t) = p2;
   v = y(t) - a2;
   k = p2 / (p2 + sigma_eps);
   a2 = a2 + k*v;
   p2 = k*sigma_eps + sigma_eta;
   vt(t) = v;
   kt(t) = k;
end

%% Plot Slide 65
x =1871:1:1970;
figure(5)
subplot(2,1,1)
plot(x(2:100),mu(2:100),'b')  % plot the time-series y in blue 'b'  
hold on
plot(x(2:100),y(2:100),'r-o')
legend('p_1 = 10^7','Observed values y_t','Interpreter','latex')

subplot(2,1,2)
plot(x(2:100),mu2(2:100),'k')  % plot the time-series y in black 
hold on
plot(x(2:100),y(2:100),'r-o')
legend('p_2 = \sigma^2_\epsilon + \sigma^2_\eta','Observed values y_t','Interpreter','latex')
%leg1 = legend('$\bar{x}$','$\tilde{x}$','$\hat{x}$');
%set(leg1,'Interpreter','latex');
%set(leg1,'FontSize',17);