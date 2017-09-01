clear all
close all
clc

load ABFe_frequency.mat
% beginning with visualization of the smoothness
figure()
plot(omega,smooth(squeeze(B(3,3,:)),0.03))
title('radiation')
figure()
plot(omega,squeeze(A(3,3,:)))
title('added mass')
figure()
plot(omega,abs(Fe(:,3)))
title('excitation')

% 
NemohData.W = omega';
NemohData.B_w = smooth(squeeze(B(3,3,:)),0.03);
NemohData.A_w = squeeze(A(3,3,:));
NemohData.Fe_w_real = abs(Fe(:,3));

%% beginning realization
figure()
plot(NemohData.W,NemohData.B_w,'bo')
hold on

dW = 0.1;
xx = 0:dW:60; % desired resolution
yy = spline(NemohData.W,NemohData.B_w,xx); % creating cubic slpline to fill in data gaps
plot(xx,yy,'r')
hold off

% compute the impulse response function using trapazoidal integration
h = 0.002;
t = [0:h:5];

Sum = zeros(1,length(t));
for jj = 1:length(t)
    for j = 1:length(xx)-1
        
        Trap = 1/2*(yy(j)+yy(j+1)) * dW;
        Trap = Trap * cos(xx(j)*t(jj));
        Sum(jj) = Sum(jj)+Trap;
    end    
end
IRF = 2/pi*Sum;

% using imp2ss to determine the impulse response function
[G] = imp2ss(IRF*h,h); % propper scaling based on Nathan Tom's discussions

% analyzing the stability of A_r
size(G.a)
lambda = eig(G.a);
MAX = max(real(lambda))
AVG = mean(real(lambda))

% looking at the comparison between the original frequency response
H = freqresp(G,xx);
for j = 1:length(xx)
    Mag(j) = abs(real(H(:,:,j)));
end

figure()
plot(xx,yy,'b',xx,Mag,'r')
legend('NEMOH','realization')
title('F.R.F.')

% looking at the impulse response of the system
[Y,time] = impulse(G,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('I.R.F.')

%% Reducing the order of the system

NumStates = 5; % Number of states for the reduced SS realization
MRTYPE = 1; % stands for model reduction type
[AM,BM,CM,DM,TOTBND,HSV] = balmr(G.a,G.b,G.c,G.d,MRTYPE,NumStates);

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

A_r = AM;
B_r = BM;
C_r = CM;
D_r = DM;
save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')

stability = eig(A_r)


% Analyzing the quality of the regression
K_r = IRF';
K_rbar = mean(IRF);
K_rtilde = Y;

Rsquared = 1-sum((K_r-K_rtilde).^2)/sum((K_r-K_rbar).^2)




