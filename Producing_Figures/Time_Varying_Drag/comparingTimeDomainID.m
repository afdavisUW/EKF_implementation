clear all
close all
clc

% Plotting comparisons of time domain drag identificaitons at varous
% frequencies

% load('EKF4.mat')
% EKF4 = EKF;


figure()
% first plot
load('EKF4_1_1Ts.mat')
EKF4.x.hat = temp;
EKF4.time.hat = tempTime;
load('EKFSS1_1Ts.mat')

subplot(3,1,1)
hold on
plot(Time.t,x.noisy(:,3),'k-.')
plot(EKF4.time.hat,EKF4.x.hat(3,:),'b')
hold off

axis([0,30,0.2,1])
ylabel('C_D')
title('Drag Period = T_p')

clear all
% second plot
load('EKF4.mat')
% EKF4.x.hat = temp;
% EKF4.time.hat = tempTime;
load('EKFSS.mat')

subplot(3,1,2)
hold on
plot(Time.t,x.noisy(:,3),'k-.')
plot(EKF.time.hat,EKF.x.hat(3,:),'b')
axis([0,30,0.2,1])
hold off
ylabel('C_D')
title('Drag Period = 3*T_p')


clear all
% third plot
load('EKF4_1_5Ts.mat')
EKF4.x.hat = temp;
EKF4.time.hat = tempTime;
load('EKFSS1_5Ts.mat')

subplot(3,1,3)
hold on
plot(Time.t,x.noisy(:,3),'k-.')
plot(EKF4.time.hat,EKF4.x.hat(3,:),'b')
hold off
axis([0,30,0.2,1])
title('Drag Period = 5*T_p')
ylabel('C_D')
xlabel('Simulation Time [s]')
legend ('True Value','4^{th} order SS')





