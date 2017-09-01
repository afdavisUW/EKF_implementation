clear all
close all
clc
format compact
set(0,'DefaultFigurePosition','factory')
% this code will be used to do independent processing on the SS realization
% of the ABFe_data.mat file

load ABFe_frequency.mat
A_BEM(1,:) = A(3,3,:);
B_BEM(1,:) = B(3,3,:);

% Specifying which data set I am working with
% dW = 0.01;
% omega = [0.1:dW:10];
omega = [[0.700357051795725,0.990454441153151,1.40071410359145,1.71551741465950,1.98090888230630,2.21472345903501,2.42610799429869,2.62049613623070,2.80142820718290,2.97136332345945,3.13209195267317,3.43103482931899,3.70594117600374,3.96181776461260,4.20214231077435,4.42944691807002,4.95227220576575,5.42494239600754,6.26418390534633,7.00357051795725,7.67202711152665,8.28673639015988,8.85889383614004,9.39627585801950,9.90454441153151;]];
W = omega;
h = 0.01;
t = [0:h:5];



% computing the impulse response function
% writing my own trapazoidal integration scheme
Sum = zeros(1,length(t));


for jj = 1:length(t)
    for j = 1:length(W)-1
                dW = W(j+1)-W(j);
        Trap = 1/2*(B_BEM(j)+B_BEM(j)) * dW;
        Trap = Trap * cos(W(j)*t(jj));
        Sum(jj) = Sum(jj)+Trap;
    end    
end
IRF = 2/pi*Sum;

% figure
% plot(t,IRF)

% using imp2ss to determine the impulse response function
[G] = imp2ss(IRF*h,h); % propper scaling based on Nathan Tom's

global A_r B_r C_r D_r
A_r = G.a;
B_r = G.b;
C_r = G.c;
D_r = G.d;
save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')

% analyzing the stability of A_r
size(A_r)
lambda = eig(A_r);
MAX = max(real(lambda))
AVG = mean(real(lambda))

% looking at the comparison between the original frequency response
H = freqresp(G,W);
for j = 1:length(W)
    Mag(j) = abs(real(H(:,:,j)));
end

figure()
plot(W,B_BEM,'b',W,Mag,'r')
legend('Nemoh','realization')
title('F.R.F. Numerical Artifacts Retained')

% looking at the impulse response of the system
[Y,time] = impulse(G,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('I.R.F. Numerical Artifacts Retained')


%% Fixing the numerical effects of the BEM code
% 
% 
% A_Nemoh(1,:) = A(3,3,:);
% B_Nemoh(1,:) = B(3,3,:);
% 
% % reformatting index 520:555
% dy = B_Nemoh(555)-B_Nemoh(520)
% dx = omega(555)-omega(520)
% M = dy/dx
% b = B_Nemoh(520)-M*omega(520)
% 
% % replacing the discontinuity
% for j = 1:length(omega(520:555))
%     B_Nemoh(519+j) = M*omega(519+j)+b;
% end
% 
% 
% % computing the impulse response function
% % writing my own trapazoidal integration scheme
% Sum = zeros(1,length(t));
% 
% 
% for jj = 1:length(t)
%     for j = 1:length(W)-1
%         
%         Trap = 1/2*(B_Nemoh(j)+B_Nemoh(j)) * dW;
%         Trap = Trap * cos(W(j)*t(jj));
%         Sum(jj) = Sum(jj)+Trap;
%     end    
% end
% IRF = 2/pi*Sum;
% 
% % figure
% % plot(t,IRF)
% 
% % using imp2ss to determine the impulse response function
% [G] = imp2ss(IRF*h,h); % propper scaling based on Nathan Tom's
% 
% global A_r B_r C_r D_r
% A_r = G.a;
% B_r = G.b;
% C_r = G.c;
% D_r = G.d;
% save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')
% 
% % analyzing the stability of A_r
% size(A_r)
% lambda = eig(A_r);
% MAX = max(real(lambda))
% AVG = mean(real(lambda))
% 
% % looking at the comparison between the original frequency response
% H = freqresp(G,W);
% for j = 1:length(W)
%     Mag(j) = abs(real(H(:,:,j)));
% end
% 
% figure()
% plot(W,B_Nemoh,'b',W,Mag,'r.')
% legend('Nemoh','realization')
% legend('Nemoh','realization')
% title('F.R.F. Numerical Artifacts Removed')
% 
% % looking at the impulse response of the system
% [Y,time] = impulse(G,t);
% figure()
% plot(time,Y,'r.',t,IRF,'b')
% title('IRF with Numerical Artifacts removed')

%% Trying to reduce the order of the state space realization

NumStates = 4; % Number of states for the reduced SS realization
MRTYPE = 1; % stands for model reduction type
[AM,BM,CM,DM,TOTBND,HSV] = balmr(A_r,B_r,C_r,D_r,MRTYPE,NumStates);

% analyzing the stability of A_r
size(AM)
lambda = eig(AM);
MAX = max(real(lambda))
AVG = mean(real(lambda))

clear H Mag
reduced_model = ss(AM,BM,CM,DM);
% looking at the comparison between the original frequency response
H = freqresp(reduced_model,W);
for j = 1:length(W)
    Mag(j) = abs(real(H(1,1,j)));
end

figure()
plot(W,B_BEM,'b',W,Mag,'r.')
legend('Nemoh','realization')
title('F.R.F. Reduced System')
ylabel('damping [kg/s]')
xlabel('frequency [rad/s]')

% looking at the impulse response of the system
[Y,time] = impulse(reduced_model,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('IRF of the reduced system')
ylabel('damping [kg/s]')
xlabel('frequency [rad/s]')
legend('Realization','Nemoh')


clear A_r B_r C_r D_r
global A_r B_r C_r D_r
A_r = AM;
RadStability = eig(A_r)
B_r = BM;
C_r = CM;
D_r = DM;
save('ABCD_BEM_sphere.mat','A_r','B_r','C_r','D_r')

