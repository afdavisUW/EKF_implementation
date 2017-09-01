clear all
close all
clc
tic

rng(1)

% define global variables 
global Param Time Func Data %Test

%% Defining Model Parameters and Datasets

% Main Function Options
READ_NEW_DATA = [0]; % ? 0->No, 1->Yes;
mooringBool = [0]; % 0-> taught & 1-> catenary
simulationBool = [0]; % 0-> Use experiemental data, 1-> simulate data
Param.RTSBool = [0]; % 0-> Do NOT implement RTS, 1-> run RTS Smoother
% Param.Qvariation = [1]; % 0-> Q not variable, 1-> Q is estimated

load ABCD_r_SSreal.mat
Param.A_r = A_r;
Param.B_r = B_r;
Param.C_r = C_r;
Param.D_r = D_r;

% parameters
samplingFreq.motion = 200;
samplingFreq.wave = 128;
samplingFreq.loadCell = 128;

% name.sensorFile = 'Run_117.csv';
% name.waveFile = 'Run_117.txt';
% name.motionFile = 'Run_117_6D.tsv'; %NOTES: this file needs to be opened
name.sensorFile = 'Run_113.csv';
name.waveFile = 'Run_113.txt';
name.motionFile = 'Run_113_6D.tsv'; %NOTES: this file needs to be opened
% name.sensorFile = 'Run_099.csv';
% name.waveFile = 'Run_099.txt';
% name.motionFile = 'Run_099_6D.tsv'; %NOTES: this file needs to be opened
Param.wavePeriod = 3; % for run 113
% Param.wavePeriod = 2.06; % for run 099
% Param.wavePeriod = 4.74; % for run 117
peakFrequency = 1/Param.wavePeriod;
waveHs = 0.2; % valid for 113 and 117
% waveHs = 0.1 % valid for 099
waveLength = 9.81*Param.wavePeriod/(2*pi);


% First defining parameters that are changed regularly
Time.dT = 1/samplingFreq.motion;
Time.simStepsPerDT = 50;
Time.dtSim = Time.dT/Time.simStepsPerDT;
Time.start = 10; % min 1 period into timespan
Time.end = 13;% max 1 period from the end of the timespan
Time.t = [Time.start:Time.dT:Time.end];

% Number of states and control inputs
Param.sizeRadSS = 5;
Param.numIDStates = 1;
Param.numStates = 2+Param.numIDStates+Param.sizeRadSS; % 2 for dynamics 
Param.numOutputs = 2; 

% Noise Parameter Values
Param.Q = 0.01^2;
Param.R = (0.005)^2*eye(Param.numOutputs);
% Param.Q = 0;
% Param.R = 0*eye(Param.numOutputs);
Param.G = [0,0,1,0,0,0,0,0]';


% process noise covariance matrix
Param.GQGt = Param.G*Param.Q*Param.G';
% Param.GQGt = Param.Q.*eye(Param.numStates);

% Initial Condtions 
Param.IC.est = [0, 0, 2.5,0,0,0,0,0];
Param.IC.real = [0, 0, 3, 0,0,0,0,0]; % parameters used when simulating 'data'
Param.IC.P = 1*eye(Param.numStates);

% Param structure data Properties of the experiment
Param.rho = 1000;
Param.g = 9.81;
Param.mass = 9.582; % 1/10 scale FPS
Param.upperPanelLength = 0.5/(sqrt(2)+1);
Param.lowerPanelLength = 0.35/(sqrt(2)+1);
Param.avgPanelLength = (Param.upperPanelLength+Param.lowerPanelLength)/2;
Param.lowerArea = 2*(1+sqrt(2))*Param.lowerPanelLength^2;
Param.upperArea = 2*(1+sqrt(2))*Param.upperPanelLength^2;
Param.avgArea = 2*(1+sqrt(2))*Param.avgPanelLength^2;
Param.buoyDepth = 0.075; % SPECIFIED AS DEPTH FROM SURFACE

Param.waveNo = peakFrequency*2*pi; % deep water assumption implicit
Param.attenFactor = exp(-Param.buoyDepth*Param.waveNo); % calculated on peak frequency
% Param.attenFactor = 1;

% parameters that form the odes
Param.addedMass = 20; % Ainf from NEMOH
Param.exciteTerm = 1800;
% Param.fluidInertia = Param.rho*Param.dispVol+Param.addedMass;
Param.weight = Param.mass*Param.g;
Param.totalInertia = Param.mass+Param.addedMass;
% Param.hydrostaticCoeff = Param.lowerArea*Param.rho*Param.g;
% Param.hydrostaticCoeff = Param.avgArea*Param.rho*Param.g;
Param.hydrostaticCoeff = Param.upperArea*Param.rho*Param.g;
Param.heaveMoorStiff = 400; % [N/m]
Param.heaveMoorStiff = 500;

% Simple Mooring Coefficients
Param.dispVol = 11322077e-9; % from solidworks
Param.Delta = Param.dispVol*Param.rho;
Param.dragFactor = 1/2*Param.rho*Param.upperArea;

% ODE solver options
% Param.options = odeset('RelTol',1e-8,'AbsTol',1e-9);
Param.options = odeset();


%% Read motion, wave, and load-cell data from dataFiles folder or simulate 
% as necessary
if simulationBool == 0
    if READ_NEW_DATA == 0
        load IncidentWaveData.mat
            Func = generateFunctions();  % lowercase functions -> nonlinear
    Func.fODE = @fNonlin;
    y.noise = 0*sqrt(Param.R)*(randn(length(Time.t), Param.numOutputs))';
    for k = 1:length(Time.t) % NOTE: noisy output is generated using noisy states (double noise)
        Data.motion.FPS.heave(k+2001,1) = Data.motion.FPS.heave(k+2001,1) + y.noise(1,k);
        Data.motion.FPS.heaveVel(k+2001,1) = Data.motion.FPS.heaveVel(k+2001,1)+y.noise(2,k);
    end 

    else
        readIncidentData(name, samplingFreq, mooringBool)
    end
%     Measurements = [Data.motion.FPS.heave];
    Measurements = [Data.motion.FPS.heave(1:end-1), Data.motion.FPS.heaveVel];
else
    
%     Param.heaveMoorStiff = 420
%     display('Placeholder Heave Mooring stiffness used!')
%     functions need to be generated early for use here
    Func = generateFunctions();  % lowercase functions -> nonlinear
    Func.fODE = @fNonlin;
    load IncidentWaveData.mat

%     dt = 1/samplingFreq.wave;
%     Data.wave.time = [Time.start-Param.wavePeriod:dt:1.5*Time.end+Param.wavePeriod];
%     Data.wave.incidentFPS = waveHs/2*sin(2*pi*1/Param.wavePeriod*Data.wave.time);

    % taking numerical derivitives
%     span = 0.01; % span percentage for smoothing
%     Data.wave.incidentFPSVel = diff(Data.wave.incidentFPS)/dt;
    % Data.wave.incidentFPSVel = smooth(Data.wave.time(1:length(tempVel)),tempVel,span,'rloess');
%     Data.wave.incidentFPSAcc = diff(Data.wave.incidentFPSVel)/dt;
    % Data.wave.incidentFPSAcc = smooth(Data.wave.time(1:length(tempAcc)),tempAcc,span,'rloess');

    Data.motion.time = Time.t;
    [~,trueData] = ode45(Func.fODE, Time.t, Param.IC.real);
    x.true = trueData;


    for k = 1:length(x.true)
        y.true(:,k) = Func.h(x.true(k,:));
    end
    x.noise = Param.G*sqrt(Param.Q)*randn(1, length(Data.motion.time));
% x.noise
%     load xtest.mat
%     x.noise = tempx(:,1:length(Time.t));
    y.noise = sqrt(Param.R)*(randn(length(Time.t), Param.numOutputs))';
%     load ytest.mat
%     y.noise = tempy(:,1:length(Time.t));
    x.noisy = x.true + x.noise';
    for k = 1:length(x.true) % NOTE: noisy output is generated using noisy states (double noise)
        y.noisy(:,k) = Func.h(x.noisy(k,:)) + y.noise(:,k);
%         y.temp(:,k) = Func.h(x.noisy(k,:));
    end 
    Measurements = y.noisy';
    Data.motion.FPS.heave = y.noisy(1,:);
    Data.motion.FPS.heaveVel = y.noisy(2,1:end-1);
end




% Uncomment to show incident wave data
    figure()    
    hold on
    plot(Data.wave.time,Data.wave.incidentFPS,'r')
    plot(Data.motion.time,Data.motion.FPS.heave,'b')    
    ylabel('displacements [m]')
    xlabel('time [S]')
    legend('Incident wave','FPS heave')
    title('Heave')
    hold off
    
    figure()
    plot(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,'r',Data.motion.time(1:end-1),Data.motion.FPS.heaveVel,'b')
    title('wave Velocity')
    legend('Incident Wave','FPS vel')
% 
%     figure()
%     plot(Data.wave.time(1:end-2),Data.wave.incidentFPSAcc)
%     title('wave acceleration')

%% Specify Function Handles
Func = generateFunctions();  % lowercase functions -> nonlinear
% Func.f = @f; % f is only used in fNonlin ODE function
% Func.h = @h;
% Func.F = @Jacob_f;
% Func.H = @Jacob_h;
Func.fODE = @fNonlin;


% Run the EKF and the RTS smoother
[val,index] = min(abs(Data.motion.time-Time.start)); % determine index for start time


%% Running simulation with estimated parameter
% Param.IC.identified = [0.02, 0, EKF.x.hat(3,end),EKF.x.hat(4,end),0.043,0,EKF.x.hat(7,end)];
[tempval,startIndex] = min(abs(Data.motion.time-Time.start));
[tempval,endIndex] = min(abs(Data.motion.time-Time.end));

% drags = [0,.35,1,2,3,4];
drags = [2:.1:4];

for j = 1:length(drags)
% Param.IC.identified = [x.true(startIndex,1), x.true(startIndex,1), drags(j),0,0,0,0,0];
Param.IC.identified = [Data.motion.FPS.heave(startIndex), Data.motion.FPS.heaveVel(startIndex), drags(j),0,0,0,0,0];
% Param.IC.identified = [0, 0, drags(j),0,0,0,0,0];

display(num2str(drags(j)))
[~,Xsimulated] = ode45(Func.fODE, Time.t, Param.IC.identified);

% for simulated data
% Error1(j) = sum(abs(Data.motion.FPS.heave(startIndex:endIndex)-Xsimulated(:,1)'));
% Error2(j) = sum(abs(Data.motion.FPS.heaveVel(startIndex:endIndex-1)-Xsimulated(1:end-1,2)'));

Error1(j) = sum(abs(Data.motion.FPS.heave(startIndex:endIndex)-Xsimulated(:,1)));
Error2(j) = sum(abs(Data.motion.FPS.heaveVel(startIndex:endIndex)-Xsimulated(1:end,2)));
Error = Error1+Error2;
end

[val,minIndex] = min(Error)
DRAG = drags(minIndex)

% Param.IC.identified = [0, 0, drags(j),0,0,0,0,0];
Param.IC.identified = [Data.motion.FPS.heave(startIndex), Data.motion.FPS.heaveVel(startIndex), drags(minIndex),0,0,0,0,0];

[~,Xsimulated] = ode45(Func.fODE, Time.t, Param.IC.identified);


for k = 1:length(Time.t)
    [U1(k),U2(k),U3(k)] = particleDynamics(Time.t(k));
end 
figure()
% subplot(2,1,1)
plot(Time.t,Xsimulated(:,1),Data.motion.time(startIndex:endIndex),Data.motion.FPS.heave(startIndex:endIndex),Time.t,U1)
title('Simulation Vs Experiment')
legend('Sim','Exper','Wave')
% xlabel('Time [s]')
ylabel('Heave [m]')


figure()
plot(Time.t,Xsimulated(:,2),Data.motion.time(startIndex:endIndex),Data.motion.FPS.heaveVel(startIndex:endIndex),Time.t,U2)
title('Sim V Exp Velocity')
legend('Sim','Exp','Wave')
ylabel('Vel m/s')
xlabel('time s')


