%%  TIME SERIES ECONOMETRICS
%
%   ASSIGNMENT 2: Simulation
%   Plotting figure 2.4 
%
%   Charlotte Taman, Femke Vedder, Rose Barzilai, Zuzana Leova (Group 1) clear all
%% 0. Clean Workspace and Command window
clear all        %clear workspace
clc              %clear command window

%% 1. Setup
load('nile.mat')

%% 2. Kalman filter

% Initialisation
a = 0;
P = 10^7;
Ht = 15099;
Qt = 1469.1;
yt = data;
T = size(yt,1);

% Calculate values of the Kalman filter
for t = 1:size(data,1)
   v = yt(t,1) - a;
   K = P / (P + Ht);
   a = a + K*v;
   P = K*Ht + Qt;
   vt(t,1) = v;
   Kt(t,1) = K;
   Pt(t,1) = P;
   at(t,1) = a;
end
K = P / (P + Ht); 


epsilon = yt - at;
eta = at(2:100,1) - at(1:99,1);

%% 3. Simulation smoothing: 
%For simulation smoothing we used paragraph 4.9 from the book of Durbin and
%Koopman.
Ft = Pt + Ht;
Lt = ones(T,1) - Kt;
r = 0;
N = 0;

for t=T:(-1):1
    Ct(t,1) = Ht - Ht*(1/Ft(t,1) + Kt(t,1)^2*N);
    Wt(t,1) = Ht*(1/Ft(t,1) - Kt(t,1)*N*Lt(t,1));
    N = 1/Ft(t,1) + Wt(t,1)^2/Ct(t,1) + Lt(t,1)^2*N;
    d = normrnd(0,sqrt(Ct(t,1)));
    con_eps(t,1) = d + Ht*(vt(t,1)/Ft(t,1) - Kt(t,1)*r);
    r = vt(t,1)/Ft(t,1) - Wt(t,1)*d/Ct(t,1) + Lt(t,1)*r;
end

con_at = yt - con_eps;
con_eta = con_at(2:100,1) - con_at(1:99,1);

%% 4. E(alpha|Y_n) 
trans_L = Lt';
r=0;
N=0;

for t = 100:-1:1
    alpha(t,1) = at(t,1) + Pt(t,1)*r;
    r = (1/Ft(t,1))*vt(t,1) + Lt(t,1)*r;
    N = (1/Ft(t,1)) + trans_L(1,t)*N*Lt(t,1);
    vt(t,1) = Pt(t,1) - Pt(t,1)*N*Pt(t,1);
end

%% 5. Simulation at and yt

S_mu = alpha(2,1);
for t = 1:100
    S_eta = normrnd(0,sqrt(Qt));
    S_eps = normrnd(0,sqrt(Ht));
    S_yt(t,1) = S_mu + S_eps;
    S_mu = S_mu + S_eta;
    S_att(t,1) = S_mu;
end

%% 6. Kalman Filter for S_yt
S_a=0;
S_P=10^7;

for t = 1:size(S_yt,1)
   S_v = S_yt(t,1) - S_a;
   S_K = S_P / (S_P + Ht);
   S_a = S_a + S_K*S_v;
   S_P = S_K*Ht + Qt;
   S_vt(t,1) = S_v;
   S_Kt(t,1) = S_K;
   S_Pt(t,1) = S_P;
   S_at(t,1) = S_a;
end

S_Ft = S_Pt + Ht;
S_Lt = ones(T,1) - S_Kt;

%% 7. Smoothing E(alpha|Y_n)
S_trans_L = S_Lt';
r = 0;
N = 0;

for t = 100:-1:1
   S_alpha(t,1) = S_at(t,1) + S_Pt(t,1)*r;
   r = (1/S_Ft(t,1))*S_vt(t,1) + S_Lt(t,1)*r;
   N = (1/S_Ft(t,1)) + S_trans_L(1,t)*N*S_Lt(t,1);
   S_vt(t,1) = S_Pt(t,1) - S_Pt(t,1)*N*S_Pt(t,1);
end 

con_alpha = alpha - S_alpha + S_att;

con_eta2 = con_alpha(2:100,1) - con_alpha(1:99,1);

%% plots
figure(1)
ax1=subplot(2,2,1);
plot(alpha);
hold on
plot(S_att,'.');
hold off
title('(i)');
axis([0 100 400 1400]);
xticks([0 50 100])
xticklabels({'1860','1920','1970'})
subplot(2,2,2);
plot(alpha);
hold on
plot(con_alpha,'.');
title('(ii)');
axis([0 100 400 1400]);
xticks([0 50 100])
xticklabels({'1860','1920','1970'})
hold off
z = zeros(T,1);
subplot(2,2,3);
plot(con_eps,'.');
axis([0 100 -400 300]);
hold on
plot(epsilon);
plot(z,':');
axis([0 100 -400 300]);
title('(iii)');
xticks([0 50 100])
xticklabels({'1860','1920','1970'})
hold off
subplot(2,2,4);
plot(con_eta2,'.');
axis([0 100 -250 250]);
xticks([0 50 100])
xticklabels({'1860','1920','1970'})
hold on
plot(eta);
plot(z,':');
title('(iv)');
axis([0 100 -250 250]);
hold off

