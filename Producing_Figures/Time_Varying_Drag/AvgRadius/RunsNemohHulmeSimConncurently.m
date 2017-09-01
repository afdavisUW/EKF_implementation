clear all
close all
clc
set(0,'DefaultFigurePosition','factory')

%  NOTE: this is all based on a spherical buoy with a radius of 1m
%% loading Hulme realization matricies
load ABCD_r_SSreal.mat
global A_H B_H C_H D_H
A_H = A_r; B_H = B_r; C_H = C_r; D_H = D_r;
clear A_r B_r C_r D_r
% loading NEMOH realization matricies
load ABCD_BEM_sphere.mat
global A_N B_N C_N D_N
A_N = A_r; B_N = B_r; C_N = C_r; D_N = D_r;
clear A_r B_r C_r D_r

% assigning useful frequencies
num_omega = 200 ;  % Number of wave frequencies !!!!!!!!
min_omega = 0.7;
max_omega = 9.9;
omega = linspace(min_omega,max_omega,num_omega);

%% Creating a wave spectrum
frequency = omega./(2*pi);
df = (frequency(2)-frequency(1));
dW = omega(2)-omega(1);

Tp = 2; % 3 second wave period (dominant)
Hs = 1/2; % signifigant wave height

% Bretschneider spectrum
for j = 1:length(frequency)
    A = frequency(j)/(Tp)^(-1);
    S(j) = 5*Hs^2/(16*(Tp)^(-1))*(1/A^5)*exp(-5/4*A^(-4));
end

% global diffOmega
% diffOmega = dW;
% creating global variables for speed
global GlobalFrequency GlobalMagnitude GlobalPhase
GlobalFrequency = frequency;
GlobalMagnitude = sqrt(2*S.*df);
PhaseTemp = rand([1,num_omega])*2*pi;
% load Phases.mat
% PhaseTemp = Phases;
GlobalPhase = PhaseTemp; % selecting the correct number of phases
% save('Phases.mat','PhaseTemp')

%% Running time domain simulations
t = [0:.01:25];

% first simulation
% defining initial conditions
NumStates = 0; % number of conventional states (not radiation)
NumStates = NumStates+1;
z = length(A_N);

q0(1,NumStates+1:NumStates+z) = zeros(1,z);

[tN,qN] = ode45('NemohRadForceODE',t,q0);
display('Nemoh')
temp1 = qN(:,end);
temp1 = abs(temp1);
mean1 = mean(temp1)
std1 = std(temp1)
max1 = max(temp1)

% second simulation
% defining initial conditions
NumStates = 0; % number of conventional states (not radiation)
NumStates = NumStates+1;
z = length(A_H);

q0(1,NumStates+1:NumStates+z) = zeros(1,z);

[tH,qH] = ode45('HulmeRadForceODE',t,q0);
display('Hulme')
temp2 = qH(:,end);
temp2 = abs(temp2);
mean2 = mean(temp2)
std2 = std(temp2)
max2 = max(temp2)

display('Percent Differences')
MeanPercent = (mean1-mean2)/mean2
StdPercent = (std1-std2)/std2
MaxPercent = (max1-max2)/max2

% visualizing the results
figure()
plot(tN,qN(:,end),tH,qH(:,end))
title('Radiation forces')
legend('Nemoh','Hulme')






