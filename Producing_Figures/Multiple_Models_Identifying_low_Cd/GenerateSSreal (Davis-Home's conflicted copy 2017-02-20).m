clear all
close all
clc

tic
% Creating SS realization of the Radiation Damping Coefficient
% omega = [0,5,10,15,20,25]; % rad/s
% b = [0,1.25,2.575,1.95,1.225,0.8250]; 
omega = [0.6,2.221,3.842,5.463,7.084,8.705,10.33,13.57,16.81,21.67,28.16];
b = [0.4927,19.36,59.1,87.86,96.86,92.78,81.89,65.78,49.03,38.61,18.67];

figure()
plot(omega,b,'bo')
hold on
dW = 0.1;
xx = 0:dW:25; % desired Resolution
yy = spline(omega,b,xx); % creating cubic spline data to fill in gaps
plot(xx,yy)
hold off

% compute the impulse response function using trapazoidal integration
h = 0.002;
t = [0:h:5];

Sum = zeros(1,length(t));
for jj = 1:length(t)
    for j = 1:length(xx)-1
        
        Trap = 1/2*(yy(j)+yy(j)) * dW;
        Trap = Trap * cos(xx(j)*t(jj));
        Sum(jj) = Sum(jj)+Trap;
    end    
end
IRF = 2/pi*Sum;
% figure
% plot(t,IRF)

% using imp2ss to determine the impulse response function
[G] = imp2ss(IRF*h,h); % propper scaling based on Nathan Tom's
% 
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
H = freqresp(G,xx);
for j = 1:length(xx)
    Mag(j) = abs(real(H(:,:,j)));
end

figure()
plot(xx,yy,'b',xx,Mag,'r')
legend('WAMIT','realization')
title('F.R.F.')

% looking at the impulse response of the system
[Y,time] = impulse(G,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('I.R.F.')

toc
%% Reducing the order

load ABCD_r_SSreal.mat

NumStates = 6; % Number of states for the reduced SS realization
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
H = freqresp(reduced_model,xx);
for j = 1:length(xx)
    Mag(j) = abs(real(H(1,1,j)));
end

figure()
plot(xx,yy,'b',xx,Mag,'r.')
legend('Nemoh','realization')
legend('Nemoh','realization')
title('F.R.F. Reduced System')

% looking at the impulse response of the system
[Y,time] = impulse(reduced_model,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('IRF of the reduced system')

clear A_r B_r C_r D_r
global A_r B_r C_r D_r
A_r = AM;
B_r = BM;
C_r = CM;
D_r = DM;
save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')

stability = eig(A_r)

toc

%% Excitation Force
omega = [0,2,4,6,8,10,12];
M = [200,192.5,160,121.25,85,61.25,41.25]; % normalized

figure
plot(omega,M,'ro')
hold on
xx = 0:0.1:12; % desired Resolution
yy = spline(omega,M,xx); % creating cubic spline data to fill in gaps
plot(xx,yy)
hold off

save('W.mat','xx')
save('Fe.mat','yy')




