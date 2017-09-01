clear all
close all
clc


load('EKFL.mat')
EKFL = EKF;
load('EKFM.mat')
EKFM = EKF;
load('EKFSS.mat')
EKFSS = EKF;


% plotting 
figure()
subplot(2,1,1)
hold on
plot(EKFSS.time.t, Measurements(:,1), 'go', EKFSS.time.hat, EKFSS.y.hat(1,:),'bx')
plot(EKFL.time.hat,EKFL.y.hat(1,:),'r+')
plot(EKFM.time.hat,EKFM.y.hat(1,:),'k')
title('EKF Estimations')
% hold on
% plot(EKF.time.hat, [EKF.y.hat(1,:)+EKF.sig3(1,:); EKF.y.hat(1,:)-EKF.sig3(1,:)],'b--')
hold off
% legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
legend('Data','SS','L','M')
ylabel('Position [m]')
% axis([Time.start,Time.end,-0.15,0.15])
subplot(2,1,2)
hold on
plot(EKFSS.time.t, Measurements(:,2), 'go', EKFSS.time.hat, EKFSS.y.hat(2,:),'bx')
plot(EKFL.time.hat,EKFL.y.hat(1,:),'r+')
% plot(EKFM.time.hat,EKFM.y.hat(1,:),'k')
% plot(Data.motion.time(1:end), Data.motion.FPS.heaveVel, 'k', EKF.time.hat, EKF.y.hat(2,:),'b')
% plot(EKF.time.hat, [EKF.y.hat(2,:)+EKF.sig3(2,:); EKF.y.hat(2,:)-EKF.sig3(2,:)],'r--')
hold off
% axis([Time.start,Time.end,-0.3,0.3])
xlabel('Time [s]')
ylabel('Velocity [m]')
% legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
legend('Data','SS','L','M')



% displayig identified drag coefficients
figure()
hold on
plot(EKFSS.time.hat,EKFSS.x.hat(3,:),'b')
plot(EKFL.time.hat,EKFL.x.hat(3,:),'r')
% plot(EKFM.time.hat,EKFM.x.hat(3,:),'k')

title('EKF Identified Drag')
hold on
plot(EKFSS.time.hat, [EKFSS.x.hat(3,:)+EKFSS.sig3(3,:); EKFSS.x.hat(3,:)-EKFSS.sig3(3,:)],'b--')
plot(EKFL.time.hat, [EKFL.x.hat(3,:)+EKFL.sig3(3,:); EKFL.x.hat(3,:)-EKFL.sig3(3,:)],'r--')
% plot(EKFM.time.hat, [EKFM.x.hat(3,:)+EKFM.sig3(3,:); EKFM.x.hat(3,:)-EKFM.sig3(3,:)],'k--')

hold off
legend('SS','L','SS 3\sigma','L 3\sigma','Location','SouthWest')
ylabel('Drag')

