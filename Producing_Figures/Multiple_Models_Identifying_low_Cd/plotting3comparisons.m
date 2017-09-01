clear all
close all
clc


load('EKFL.mat')
EKFL = EKF;
load('EKFM.mat')
EKFM = EKF;
load('EKF1.mat')
EKF1 = EKF;
load('EKF2.mat')
EKF2 = EKF;
load('EKF3.mat')
EKF3 = EKF;
load('EKF4.mat')
EKF4 = EKF;
load('EKFSS.mat')
EKFSS = EKF;


% plotting 
figure()
subplot(2,1,1)
hold on
plot(EKFSS.time.t, Measurements(:,1), EKF4.time.hat,EKF4.y.hat(1,:), EKF3.time.hat,EKF3.y.hat(1,:), EKF2.time.hat,EKF2.y.hat(1,:), EKF1.time.hat,EKF1.y.hat(1,:), EKFL.time.hat,EKFL.y.hat(1,:), EKFM.time.hat,EKFM.y.hat(1,:),'k.')
title('EKF Estimations')
% hold on
% plot(EKF.time.hat, [EKF.y.hat(1,:)+EKF.sig3(1,:); EKF.y.hat(1,:)-EKF.sig3(1,:)],'b--')
hold off
% legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
legend('Data','SS4','SS3','SS2','SS1','L','M')
ylabel('Position [m]')
% axis([Time.start,Time.end,-0.15,0.15])
subplot(2,1,2)
hold on
plot(EKFSS.time.t, Measurements(:,2), 'g')
plot(EKFSS.time.hat, EKFSS.y.hat(2,:),'b')
plot(EKF2.time.hat,EKF2.y.hat(2,:),'m')
plot(EKF1.time.hat,EKF1.y.hat(2,:),'c')
plot(EKFL.time.hat,EKFL.y.hat(2,:),'r')
plot(EKFM.time.hat,EKFM.y.hat(2,:),'k')
% plot(Data.motion.time(1:end), Data.motion.FPS.heaveVel, 'k', EKF.time.hat, EKF.y.hat(2,:),'b')
% plot(EKF.time.hat, [EKF.y.hat(2,:)+EKF.sig3(2,:); EKF.y.hat(2,:)-EKF.sig3(2,:)],'r--')
hold off
% axis([Time.start,Time.end,-0.3,0.3])
xlabel('Time [s]')
ylabel('Velocity [m]')
% legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
legend('Data','SS','SS2','SS1','L','M')



% displayig identified drag coefficients
figure()
subplot(2,1,1)
hold on
% plot(EKFSS.time.hat,EKFSS.x.hat(3,:),'b')
% plot(EKFL.time.hat,EKFL.x.hat(3,:),'r')
% plot(EKF2.time.hat,EKF2.x.hat(3,:),'g')
% plot(EKF1.time.hat,EKF1.x.hat(3,:),'c')
plot(EKF4.time.hat,EKF4.x.hat(3,:), EKF3.time.hat,EKF3.x.hat(3,:), EKF2.time.hat,EKF2.x.hat(3,:),  EKFL.time.hat,EKFL.x.hat(3,:), EKF1.time.hat,EKF1.x.hat(3,:), EKFM.time.hat,EKFM.x.hat(3,:))
plot(EKFSS.time.hat,(0.35)*ones(1,length(EKFSS.time.hat)),'k-.')
% plot(EKFSS.time.hat,(0.35+0.1*0.35)*ones(1,length(EKFSS.time.hat)),'k--')
% plot(EKFSS.time.hat,(0.35-0.1*0.35)*ones(1,length(EKFSS.time.hat)),'k--')


% plot(EKFM.time.hat,EKFM.x.hat(3,:),'k')

title('Identified Drag')
legend('SS4','SS3','SS2','Linear','SS1','Morrison','location','NorthWest','orientation','Horizontal')
axis([0,10,-0.1,2.5])
ylabel('C_D')
% hold on
% plot(EKFSS.time.hat, [EKFSS.x.hat(3,:)+EKFSS.sig3(3,:); EKFSS.x.hat(3,:)-EKFSS.sig3(3,:)],'b--')
% plot(EKFL.time.hat, [EKFL.x.hat(3,:)+EKFL.sig3(3,:); EKFL.x.hat(3,:)-EKFL.sig3(3,:)],'r--')
% plot(EKF2.time.hat, [EKF2.x.hat(3,:)+EKF2.sig3(3,:); EKF2.x.hat(3,:)-EKF2.sig3(3,:)],'g--')
% plot(EKF1.time.hat, [EKF1.x.hat(3,:)+EKF1.sig3(3,:); EKF1.x.hat(3,:)-EKF1.sig3(3,:)],'c--')
% plot(EKFM.time.hat, [EKFM.x.hat(3,:)+EKFM.sig3(3,:); EKFM.x.hat(3,:)-EKFM.sig3(3,:)],'k--')
hold off
% legend('SS','L','SS 3\sigma','L 3\sigma','Location','SouthWest')
ylabel('Drag Coefficient')


subplot(2,1,2)
hold on
plot(EKF4.time.hat,EKF4.x.hat(3,:), EKF3.time.hat,EKF3.x.hat(3,:), EKF2.time.hat,EKF2.x.hat(3,:), EKFL.time.hat,EKFL.x.hat(3,:))
plot(EKFSS.time.hat,(0.35)*ones(1,length(EKFSS.time.hat)),'k-.')
plot(EKFSS.time.hat,(0.35+0.1*0.35)*ones(1,length(EKFSS.time.hat)),'k--')
plot(EKFSS.time.hat,(0.35-0.1*0.35)*ones(1,length(EKFSS.time.hat)),'k--')
hold off
axis([0,10,0.1,0.5])
legend('SS4','SS3','SS2','Linear','location','SouthWest')
title('Magnification of selected models')
xlabel('Simulation Time [s]')
ylabel('Drag Coefficient')

