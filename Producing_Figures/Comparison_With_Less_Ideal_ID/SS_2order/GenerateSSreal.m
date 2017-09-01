clear all
close all
clc

tic
% Creating SS realization of the Radiation Damping Coefficient
% omega = [0,5,10,15,20,25]; % rad/s
% b = [0,1.25,2.575,1.95,1.225,0.8250]; 
% omega = [0.6,2.221,3.842,5.463,7.084,8.705,10.33,13.57,16.81,21.67,28.16];
% b = [0.4927,19.36,59.1,87.86,96.86,92.78,81.89,65.78,49.03,38.61,18.67];
omega = [0.600000000000000,1.88333333333333,3.16666666666667,4.45000000000000,5.73333333333333,7.01666666666667,8.30000000000000,9.58333333333333,10.8666666666667,12.1500000000000,13.4333333333333,14.7166666666667,16.0000000000000,17.2833333333333,18.5666666666667,19.8500000000000,21.1333333333333,22.4166666666667,23.7000000000000,24.9833333333333,26.2666666666667,27.5500000000000,28.8333333333333,30.1166666666667,31.4000000000000];

b = [0.492670700000000,12.7727500000000,42.2720300000000,72.2761800000000,90.6180000000000,96.7930000000000,94.6106100000000,87.4484600000000,77.2938000000000,53.5029200000000,67.5231200000000,53.7529800000000,58.7994100000000,46.3644600000000,28.3712200000000,37.1208500000000,31.3246700000000,30.0966900000000,26.5735400000000,24.4216500000000,20.9871300000000,19.2334000000000,18.0881900000000,14.6836200000000,13.7246400000000];

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
% omega = [0,2,4,6,8,10,12];
% omega = [0.600000000000000,1.88333333333333,3.16666666666667,4.45000000000000,5.73333333333333,7.01666666666667,8.30000000000000,9.58333333333333,10.8666666666667,12.1500000000000,13.4333333333333,14.7166666666667,16.0000000000000,17.2833333333333,18.5666666666667,19.8500000000000,21.1333333333333,22.4166666666667,23.7000000000000,24.9833333333333,26.2666666666667,27.5500000000000,28.8333333333333,30.1166666666667,31.4000000000000];
M = [2066.23200000000;1891.86900000000;1578.74900000000;1239.56400000000;949.271200000000;724.490500000000;556.803000000000;431.938800000000;337.947100000000;258.575700000000;218.990300000000;172.564700000000;159.828400000000;126.635900000000;92.8322500000000;91.1462900000000;75.4302200000000;66.6539500000000;57.5478600000000;48.4123500000000;41.7192200000000;35.8349800000000;31.4644300000000;24.8084200000000;22.4794500000000];
% M = [200,192.5,160,121.25,85,61.25,41.25]; % normalized

figure
plot(omega,M,'ro')
hold on
xx = 0:0.1:31; % desired Resolution
yy = spline(omega,M,xx); % creating cubic spline data to fill in gaps
plot(xx,yy)
hold off

save('W.mat','xx')
save('Fe.mat','yy')




