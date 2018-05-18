%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 6: LLM and Nile data
%
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1)
%   March 2018 

%% 0. Clean Workspace and Command Window

clear all        %clear workspace
clc              %clear command window

%% 1. Read Data
load('nile.mat')

%% 2. Parameter Initialisation
    
    sigma_eps = 15099;
    sigma_eta = 1469;
    yt = data;
    T = size(yt,1);
    
%% 3. Filtered state and filtered state variance SIS N=100

N=100;
mu = yt(1,1)*ones(1,N);
p = zeros(T,N);
w(1,:) = ones(1,N)/N;

for t=1:T
for n=1:N
    mu(t+1,n) = mu(t,n) + normrnd(0,sqrt(sigma_eta)); %state equation, alpha=mu in LLM
    p(t,n) = 1/sqrt(2*pi*sigma_eps) * exp(-(yt(t,1) - mu(t,n))^2/(2*sigma_eps)); %density p(y_t|mu_t) 
end

w(t+1,:) = w(t,:).*p(t,:); %corresponsing weights
w_norm(t,:) = w(t+1,:)./sum(w(t+1,:)); %normalized weights
 
estimate_x(t,1) = sum(w_norm(t,:).*mu(t,:)); %estimates
var_estimate_x(t,1) = sum(w_norm(t,:).*((mu(t,:)-estimate_x(t,1)).^2)); %variance estimates

ESS(t,1) = 1/(sum(w_norm(t,:).*w_norm(t,:))); %effective sample size
end
      
estimate_x
%% 4. Apply Kalman Filter
% initialization
a = 0;
p = 10^7;

% kalman filter for LLM
for t = 1:T
   mu_t(t,1) = a;  
   p_t(t,1) = p;
   v = yt(t,1) - a;
   k = p / (p + sigma_eps);
   a = a + k*v;
   p = k*sigma_eps + sigma_eta;
   v_t(t,1) = v;
   K_t(t,1) = k;
end

Ft = p_t + sigma_eps;

%% 5. Bootstrap filter 

alpha = yt(1,1)*ones(1,N);
  
for t = 1:T 
for n = 1:N
    p(t,n) = 1/sqrt(2*pi*sigma_eps) * exp(-(yt(t,1) - alpha(t,n))^2/(2*sigma_eps));
end  
    w(t,:) = p(t,:)./sum(p(t,:));
for n = 1:N
    cum = cumsum([0; w(t,:).'/sum(w(t,:))]); % adds up to 1
    cum(end) = 1e3*eps + cum(end);
    [j j] = histc(rand,cum(:));
    alpha_new(t,n) = alpha(t,j);
end

alpha(t,:) = alpha_new(t,:);
ESS_bootstrap(t,1) = 1/sum((ones(1,N)/N).^2);

x(t,1) = sum(alpha(t,:))/N;
var_x(t,1) = sum((alpha(t,:)-x(t,1)*ones(1,N)).^2)/N; %variance

for n = 1:N
     alpha(t+1,n) = alpha(t,n) + normrnd(0,sqrt(sigma_eta));
end    
end


%% 6. Plots
% Nile data, Kalman filter and SIS (filtered state)
figure(1);
t = (1871:1970)';
plot(datenum(t,1,1),yt);
hold on
plot(datenum(t,1,1),[mu_t(2:100);a]);
plot(datenum(t,1,1),estimate_x);
legend('Nile data','Kalman filter','SIS')
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 400 1500]);
hold off

Variance of the kalman filter and SIS with N=100 (filtered state)
figure(2);
t = (1871:1970)';
plot(datenum(t,1,1),p_t);
hold on
plot(datenum(t,1,1),var_estimate_x);
legend('Kalman filter','SIS')
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 0 30000]);
hold off

% ESS
figure(3)
subplot(2,1,1)
t = (1871:1970)';
plot(datenum(t,1,1),ESS);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 0 110]);
hold off

subplot(2,1,2)
t = (1871:1970)';
plot(datenum(t,1,1),ESS_bootstrap);
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 0 110]);
hold off

% Nile data, kalman filter and bootstrap filter
figure(4);
t = (1871:1970)';
plot(datenum(t,1,1),yt);
hold on
plot(datenum(t,1,1),[mu_t(2:100);a]);
plot(datenum(t,1,1),x);
legend('Nile data','Kalman filter','Bootstrap filter')
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 400 1500]);
hold off

% Variance Kalman filter and bootstrap filter
figure(5);
t = (1871:1970)';
plot(datenum(t,1,1),p_t);
hold on
plot(datenum(t,1,1),var_x);
legend('Kalman filter','Bootstrap filter')
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1970,1,1) 0 30000]);
hold off

% Bootstrap filter with N=15 particle paths over time from t=1 upto
% t=5,10,15,20
figure(6);
subplot(2,2,1);
t = (1871:1890)';
plot(datenum(t,1,1),mu_t(2:21));
hold on
plot(datenum(t,1,1),[alpha(1:5,1:15);NaN*ones(15,15)]);
legend('Kalman filter','Bootstrap filter');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1891,1,1) 700 1500]);
hold off

subplot(2,2,2);
t = (1871:1890)';
plot(datenum(t,1,1),mu_t(2:21));
hold on
plot(datenum(t,1,1),[alpha(1:10,1:15);NaN*ones(10,15)]);
legend('Kalman filter','Bootstrap filter');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1891,1,1) 700 1500]);
hold off

subplot(2,2,3);
t = (1871:1890)';
plot(datenum(t,1,1),mu_t(2:21));
hold on
plot(datenum(t,1,1),[alpha(1:15,1:15);NaN*ones(5,15)]);
legend('Kalman filter','Bootstrap filter');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1891,1,1) 700 1500]);
hold off

subplot(2,2,4);
t = (1871:1890)';
plot(datenum(t,1,1),mu_t(2:21));
hold on
plot(datenum(t,1,1),alpha(1:20,1:15));
legend('Kalman filter','Bootstrap filter');
dateFormat = 10;
datetick('x',dateFormat);
axis([datenum(1870,1,1) datenum(1891,1,1) 700 1500]);
hold off

