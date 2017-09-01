lear all
close all
clc
tic

rng(1)

%S define global variables and frequency information
global Param Time Func Data %Test
load ABFe_frequency.mat
% load Meas.mat
global NemohData SpectrumData
NemohData.W = omega';
NemohData.f = (omega./(2*pi))';
% NemohData.B_w = smooth(squeeze(B(3,3,:)),0.03);
% NemohData.A_w = squeeze(A(3,3,:));
NemohData.Fe_w_real = abs(Fe(:,3));


%% Defining simulated spectrum data
% using a Bretscneider spectrum

Param.wavePeriod = 2.69; % for run 125
Param.waveHs = 0.3;

num = 1000; % number of sine waves summed to generate irregular waves
minOmega = 0.15;
maxOmega = 30;
frequency = linspace(minOmega/(2*pi),maxOmega/(2*pi),num);
dW = (frequency(2)-frequency(1))*2*pi;

for j = 1:length(frequency) % loop to calculate specta
    A = frequency(j)/(Param.wavePeriod^(-1));
    S(j) = 5*(Param.waveHs/2)^2/(16*(Param.wavePeriod)^(-1))*(1/A^5)*exp(-5/4*A^(-4));
%     S(j) = 0 % option preserved to generate no wave
end

SpectrumData.dW = dW;
SpectrumData.frequency = frequency;
SpectrumData.magnitude = sqrt(2*S.*dW);
SpectrumData.phase = rand([1,num])*2*pi;
SpectrumData.W = frequency.*2*pi;

% formatting frequency excitation to useful form
SpectrumData.Fe = spline(NemohData.W,NemohData.Fe_w_real,SpectrumData.W);

%% Defining Model Parameters and Datasets

% Main Function Options
READ_NEW_DATA = [0]; % ? 0->No, 1->Yes;
mooringBool = [0]; % 0-> taught & 1-> catenary
simulationBool = [1]; % 0-> Use experiemental data, 1-> simulate data
Param.RTSBool = [0]; % 0-> Do NOT implement RTS, 1-> run RTS Smoother
% Param.Qvariation = [1]; % 0-> Q not variable, 1-> Q is estimated

load ABCD_r_SSreal.mat
Param.A_r = A_r;
Param.B_r = B_r;
Param.C_r = C_r;
Param.D_r = D_r;

% sys = ss(A_r,B_r,C_r,D_r);
% t = 0:0.01:10;
% u = sin(pi*t);
% lsim(sys,u,t)
% bode(sys)
% return

% parameters
samplingFreq.motion = 200;
samplingFreq.wave = 128;
samplingFreq.loadCell = 128;

% name.sensorFile = 'Run_117.csv';
% name.waveFile = 'Run_117.txt';
% name.motionFile = 'Run_117_6D.tsv'; %NOTES: this file needs to be opened
% name.sensorFile = 'Run_113.csv';
% name.waveFile = 'Run_113.txt';
% name.motionFile = 'Run_113_6D.tsv'; %NOTES: this file needs to be opened
name.sensorFile = 'Run_125.csv';
name.waveFile = 'Run_125.txt';
name.motionFile = 'Run_125_6D.tsv'; %NOTES: this file needs to be opened
% name.sensorFile = 'Run_099.csv';
% name.waveFile = 'Run_099.txt';
% name.motionFile = 'Run_099_6D.tsv'; %NOTES: this file needs to be opened
% Param.wavePeriod = 3; % for run 113
% Param.wavePeriod = 2.06; % for run 099
% Param.wavePeriod = 4.74; % for run 117
peakFrequency = 1/Param.wavePeriod;
% waveHs = 0.2; % valid for 113 and 117
% waveHs = 0.1 % valid for 099
waveLength = 9.81*Param.wavePeriod/(2*pi);


% First defining parameters that are changed regularly
Time.dT = 1/samplingFreq.motion;
Time.simStepsPerDT = 50;
Time.dtSim = Time.dT/Time.simStepsPerDT;
Time.start = 0; % min 1 period into timespan
% Time.end = 12;
Time.end = 10;% max 1 period from the end of the timespan
Time.t = [Time.start:Time.dT:Time.end];

% determining the input wave variance
    for k = 1:length(Time.t)
        [pos,vel,acc] = particleDynamics_spectrum(Time.t(k));
        incidentWaveTemp(k) = pos;
    end
stdDevInputPos = std(incidentWaveTemp);

% Number of states and control inputs
Param.sizeRadSS = 10;
Param.numIDStates = 1;
Param.numStates = 2+Param.numIDStates+Param.sizeRadSS; % 2 for dynamics 
Param.numOutputs = 2; 

% Smoother information for linear radiation identification
Param.numStatesLin = 3;
Param.LinIC.est = [0.1,0.1,1];
Param.LinIC.P = 0.1*eye(3);


Param.R = (0.0006075)^2*eye(Param.numOutputs); % variance of the residuals
Param.GQGt = zeros(Param.numStates);
% Param.GQGt(1,1) = 0.0701^2; % variance of the wave input
% Param.GQGt(2,2) = 0.0520^2; % variance of the wave input velocity
% Param.GQGt(1,1) = 0.0666^2; % variance of the wave input
Param.GQGt(1,1) = stdDevInputPos
% Param.GQGt(2,2) = 0.0739^2; % variance of the wave input velocity
TF_mag = 5/2;
% Param.GQGt(2,2) = TF_mag*0.0666^2; % variance of the wave input velocity
Param.GQGt(2,2) = TF_mag*Param.GQGt(1,1);

% Param.GQGt(3,3) = 0.01; % variance of the drag coefficient
% Param.GQGt(3,3) = 0.001;
Param.GQGt(3,3) = 0.05
Param.GQGt(3,3) = 0.1

% process noise covariance matrix
% Param.GQGt = Param.G*Param.Q*Param.G';
% Param.GQGt = Param.Q.*eye(Param.numStates);

% Initial Condtions 
Param.IC.est = [0.1, 0.1, 1,0,0,0,0,0,0,0,0,0,0];
Param.IC.real = [0, 0, 0.35, 0,0,0,0,0,0,0,0,0,0]; % parameters used when simulating 'data'
% Param.IC.P = 0.1*eye(Param.numStates);
Param.IC.P = zeros(Param.numStates);
Param.IC.P(1,1) = 0.1; Param.IC.P(2,2) = 0.1; Param.IC.P(3,3) = 0.1;

% Param structure data Properties of the experiment
Param.linRad = 20;
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

% Param.waveNo = peakFrequency*2*pi; % deep water assumption implicit
% Param.attenFactor = exp(-Param.buoyDepth*Param.waveNo); % calculated on peak frequency
% Param.attenFactor = 1;

% parameters that form the odes
% Param.addedMass = 20; % Ainf from NEMOH
Param.addedMass = 21.15;
Param.exciteTerm = 1800;
% Param.fluidInertia = Param.rho*Param.dispVol+Param.addedMass;
Param.weight = Param.mass*Param.g;
Param.totalInertia = Param.mass+Param.addedMass;
Param.hydrostaticCoeff = Param.avgArea*Param.rho*Param.g;
% Param.hydrostaticCoeff = Param.upperArea*Param.rho*Param.g

% Param.heaveMoorStiff = 400; % [N/m]
Param.heaveMoorStiff = 500;
% Param.surgeMoorStiff = 50; % [N/m]

% Simple Mooring Coefficients
Param.dispVol = 11322077e-9; % from solidworks
Param.Delta = Param.dispVol*Param.rho;
Param.dragFactor = 1/2*Param.rho*Param.upperArea;

% ODE solver options
% Param.options = odeset('RelTol',1e-5,'AbsTol',1e-6);
Param.options = odeset();


%% Read motion, wave, and load-cell data from dataFiles folder or simulate 
% as necessary
if simulationBool == 0
    if READ_NEW_DATA == 0
        load IncidentWaveData.mat
        Func = generateFunctions();  % lowercase functions -> nonlinear
        Func.fODE = @fNonlin;
        y.noise = sqrt(Param.R)*(randn(length(Time.t), Param.numOutputs))';
%     for k = 1:length(Time.t) % NOTE: noisy output is generated using noisy states (double noise)
%         Data.motion.FPS.heave(k+2001,1) = Data.motion.FPS.heave(k+2001,1) + y.noise(1,k);
%         Data.motion.FPS.heaveVel(k+2001,1) = Data.motion.FPS.heaveVel(k+2001,1)+y.noise(2,k);
%     end 

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
%     load IncidentWaveData.mat

%     dt = 1/samplingFreq.wave;
%     Data.wave.time = [Time.start-Param.wavePeriod:dt:1.5*Time.end+Param.wavePeriod];
%     for k = 1:length(Data.wave.time)
%         [pos,vel,acc] = particleDynamics_spectrum(Data.wave.time(k));
%         Data.wave.incidentFPS(k) = pos;
%         Data.wave.incidentFPSVel(k) = vel;
%         Data.wave.incidentFPSAcc(k) = acc; 
%     end
    
%     stdDevInputPos = std(Data.wave.incidentFPS)
%     stdDevInputVel = std(Data.wave.incidentFPSVel)
%     Data.wave.incidentFPS = waveHs/2*sin(2*pi*1/Param.wavePeriod*Data.wave.time);

    % taking numerical derivitives
%     span = 0.01; % span percentage for smoothing
%     Data.wave.incidentFPSVel = diff(Data.wave.incidentFPS)/dt;
    % Data.wave.incidentFPSVel = smooth(Data.wave.time(1:length(tempVel)),tempVel,span,'rloess');
%     Data.wave.incidentFPSAcc = diff(Data.wave.incidentFPSVel)/dt;
    % Data.wave.incidentFPSAcc = smooth(Data.wave.time(1:length(tempAcc)),tempAcc,span,'rloess');

    
%% Generating numerical data

    Data.motion.time = Time.t;
    [~,trueData] = ode45(Func.fODE, Time.t, Param.IC.real);

    x.true = trueData;


    for k = 1:length(x.true)
        y.true(:,k) = Func.h(x.true(k,:));
    end
%     x.noise = Param.G*sqrt(Param.Q)*randn(1, length(Data.motion.time));
x.noise = 0
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
%     Data.motion.FPS.heaveVel = y.noisy(2,1:end);
end

return

% Uncomment to show incident wave data
%     figure()    
%     hold on
%     plot(Data.wave.time,Data.wave.incidentFPS,'r')
%     plot(Data.motion.time,Data.motion.FPS.heave,'b')    
%     ylabel('displacements [m]')
%     xlabel('time [S]')
%     legend('Incident wave','FPS heave')
%     title('Heave')
%     hold off
%     
%     figure()
%     plot(Data.wave.time(1:end),Data.wave.incidentFPSVel,'r',Data.motion.time(1:end-1),Data.motion.FPS.heaveVel,'b')
% %     plot(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,'r',Data.motion.time(1:end),Data.motion.FPS.heaveVel,'b')
%     title('wave Velocity')
%     legend('Incident Wave','FPS vel')
% % 
% %     figure()
% %     plot(Data.wave.time(1:end-2),Data.wave.incidentFPSAcc)
% %     title('wave acceleration')



%% Specify Function Handles
% Func = generateFunctions();  % lowercase functions -> nonlinear
% Func.f = @f; % f is only used in fNonlin ODE function
% Func.h = @h;
% Func.F = @Jacob_f;
% Func.H = @Jacob_h;
Func.fODE = @fNonlin;

toc
%% Run the EKF and the RTS smoother
% [val,index] = min(abs(Data.motion.time-Time.start)); % determine index for start time
% [sim] = smoother_RTS(Func, Measurements(index:end,:));
[sim] = smoother_RTS(Func, Measurements(:,:));



% storing EKF and RTS structures for plotting
EKF = sim.EKF; % excecute if RTS is active
if Param.RTSBool == 1
    RTS = sim.RTS;
end

save('EKFSS.mat','Measurements')
%%
% % Running smoother for the linear model
% FuncLin = generateFunctionsLin();
% FuncLin.f = @f_lin;
% FuncLin.h = @h_lin;
% FuncLin.F = @Jacob_f_lin;
% FuncLin.H = @Jacob_h_lin;
% FuncLin.fODE = @fNonlinLin;
% [simLin] = smoother_RTS_Lin(FuncLin, Measurements(index:end,:));
% 


%% Plot the EKF and RTS results

% First the EKF filter on its own
% Displaying the confidences intervals

% displaying heave and surge displacemetns
figure()
subplot(2,1,1)
plot(Data.motion.time, Data.motion.FPS.heave, 'k', EKF.time.hat, EKF.y.hat(1,:),'b')
title('EKF Estimations')
hold on
plot(EKF.time.hat, [EKF.y.hat(1,:)+EKF.sig3(1,:); EKF.y.hat(1,:)-EKF.sig3(1,:)],'r--')
hold off
legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
ylabel('Position [m]')
axis([Time.start,Time.end,-0.15,0.15])
subplot(2,1,2)
plot(Data.motion.time(1:end-1), Data.motion.FPS.heaveVel, 'k', EKF.time.hat, EKF.y.hat(2,:),'b')
% plot(Data.motion.time(1:end), Data.motion.FPS.heaveVel, 'k', EKF.time.hat, EKF.y.hat(2,:),'b')
hold on
plot(EKF.time.hat, [EKF.y.hat(2,:)+EKF.sig3(2,:); EKF.y.hat(2,:)-EKF.sig3(2,:)],'r--')
hold off
axis([Time.start,Time.end,-0.3,0.3])
xlabel('Time [s]')
ylabel('Velocity [m]')
legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')


% displayig identified drag coefficients
figure()
plot(EKF.time.hat,EKF.x.hat(3,:),'b')
title('EKF Identified Drag')
hold on
plot(EKF.time.hat, [EKF.x.hat(3,:)+EKF.sig3(3,:); EKF.x.hat(3,:)-EKF.sig3(3,:)],'r--')
hold off
legend('C_D','3\sigma Error Bounds')
ylabel('Drag')

% plotting the state response of the radiation realization
figure()
plot(EKF.time.hat,EKF.x.hat(4,:),EKF.time.hat,EKF.x.hat(5,:),EKF.time.hat,EKF.x.hat(6,:),EKF.time.hat,EKF.x.hat(7,:))
title('State Dynamics of radiation realization')
legend('State 1','State 2','State 3','State 4')
xlabel('Time [s]')

% RTS smoother results
if Param.RTSBool == 1 % excecute if RTS is active
%     RTS displacements
    figure()
%     subplot(2,1,1)
    plot(Data.motion.time, Data.motion.FPS.heave, 'k', RTS.time.hat, RTS.y.hat(1,:),'b')
    title('RTS Extimations')
    hold on
    plot(RTS.time.hat, [RTS.y.hat(1,:)+RTS.sig3(1,:); RTS.y.hat(1,:)-RTS.sig3(1,:)],'r--')
    hold off
    legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')
    ylabel('Heave [m]')
%     subplot(2,1,2)
%     plot(Data.motion.time, Data.motion.FPS.surge, 'k', RTS.time.hat, RTS.y.hat(2,:),'b')
%     hold on
%     plot(RTS.time.hat, [RTS.y.hat(2,:)+RTS.sig3(5,:); RTS.y.hat(2,:)-RTS.sig3(5,:)],'r--')
%     hold off
%     xlabel('Time [s]')
%     ylabel('Surge [m]')
%     legend('FPS Measured','FPS Estimated','3\sigma Error Bounds')

    % displayig identified drag coefficients
    figure()
%     subplot(2,1,1)
    plot(RTS.time.hat,RTS.x.hat(3,:),'b')
    title('RTS Identified Drag')
    hold on
    plot(RTS.time.hat, [RTS.x.hat(3,:)+RTS.sig3(3,:); RTS.x.hat(3,:)-RTS.sig3(3,:)],'r--')
    hold off
    legend('C_D','3\sigma Error Bounds')
    ylabel('Heave')
%     subplot(2,1,2)
%     xlabel('Time [s]')
%     plot(RTS.time.hat,RTS.x.hat(7,:),'b')
%     hold on
%     plot(RTS.time.hat, [RTS.x.hat(7,:)+RTS.sig3(7,:); RTS.x.hat(7,:)-RTS.sig3(7,:)],'r--')
%     hold off
%     ylabel('Surge')
%     legend('C_D','3\sigma Error Bounds')

end 


%% Running simulation with estimated parameter
% Param.IC.identified = [0.02, 0, EKF.x.hat(3,end),EKF.x.hat(4,end),0.043,0,EKF.x.hat(7,end)];
[tempval,startIndex] = min(abs(Data.motion.time-Time.start));
[tempval,endIndex] = min(abs(Data.motion.time-Time.end));


Param.IC.identified = [Data.motion.FPS.heave(startIndex), Data.motion.FPS.heaveVel(startIndex), EKF.x.hat(3,end),0,0,0,0];

if Param.RTSBool ==1
    Param.IC.identified = [0.02, -0.12, RTS.x.hat(3,end),RTS.x.hat(4,end)];
end

[tempval,startIndex] = min(abs(Data.motion.time-Time.start));
[tempval,endIndex] = min(abs(Data.motion.time-Time.end));

[~,Xsimulated] = ode45(Func.fODE, Time.t, Param.IC.identified);
for k = 1:length(Time.t)
    [U1(k),U2(k),U3(k)] = particleDynamics_spectrum(Time.t(k));
end 
figure()
% subplot(2,1,1)
plot(Time.t,Xsimulated(:,1),Data.motion.time(startIndex:endIndex),Data.motion.FPS.heave(startIndex:endIndex),Time.t,U1,Time.t,U2,Time.t,U3)
title('Simulation Vs Experiment')
legend('Sim','Exper','wave','waveVel','WaveAcc')
% xlabel('Time [s]')
ylabel('Heave [m]')


%% Modeled Forces acting on the float

noisyHeave = Data.motion.FPS.heave;
noisyVel = Data.motion.FPS.heaveVel;

filteredHeave = EKF.y.hat(1,1:50:end);
filteredVel = EKF.y.hat(2,1:50:end);
filteredX4 = EKF.x.hat(4,1:50:end);
filteredX5 = EKF.x.hat(5,1:50:end);
filteredX6 = EKF.x.hat(6,1:50:end);
filteredX7 = EKF.x.hat(7,1:50:end);

for k = 1:length(Time.t)-1
%     smooth simulation forces
%     excitationForce(k) = Param.exciteTerm*U1(k);
    excitationForce(k) = excitationForce_spec(Time.t(k));
    mooringForce(k) = -Param.heaveMoorStiff*Xsimulated(k,1);
    hydrostaticForce(k) = -Param.hydrostaticCoeff*Xsimulated(k,1);
    dragForce(k) = -Param.dragFactor*EKF.x.hat(3,end)*(Xsimulated(k,2)-U2(k))*abs(Xsimulated(k,2)-U2(k));
    radForce(k) = -(Param.C_r*[Xsimulated(k,4);Xsimulated(k,5);Xsimulated(k,6);Xsimulated(k,7)]+Param.D_r*Xsimulated(k,2));

%     noisy unfiltered data
%     excitationForceNoisy(k) = Param.exciteTerm*U1(k);
    excitationForceNoisy(k) = excitationForce_spec(Time.t(k));
    mooringForceNoisy(k) = -Param.heaveMoorStiff*noisyHeave(k);
    hydrostaticForceNoisy(k) = -Param.hydrostaticCoeff*noisyHeave(k);
    dragForceNoisy(k) = -Param.dragFactor*EKF.x.hat(3,end)*(noisyHeave(k)-U2(k))*abs(noisyHeave(k)-U2(k));
%     radForceNoisy(k) = -(Param.C_r*[x.noisy(k,4);x.noisy(k,5);x.noisy(k,6);x.noisy(k,7);x.noisy(k,8)]+Param.D_r*noisyVel(k));
    
%     filtered data
%     excitationForceFilter(k) = Param.exciteTerm*U1(k);
    excitationForceFilter(k) = excitationForce_spec(Time.t(k));
    mooringForceFilter(k) = -Param.heaveMoorStiff*filteredHeave(k);
    hydrostaticForceFilter(k) = -Param.hydrostaticCoeff*filteredHeave(k);
    dragForceFilter(k) = -Param.dragFactor*EKF.x.hat(3,end)*(filteredHeave(k)-U2(k))*abs(filteredHeave(k)-U2(k));
%     radForceFilter(k) = -(Param.C_r*[filteredX4(k);filteredX5(k);filteredX6(k);filteredX7(k);filteredX8(k)]+Param.D_r*filteredVel(k));
    
    
end

% noisy radiation force needs to becomputed seperately based on noisy
% heaveVelocity. 
sys = ss(Param.A_r,Param.B_r,Param.C_r,Param.D_r);
radForceNoisy = -lsim(sys,noisyVel,Time.t(1:end-1));
radForceFilter = -lsim(sys,filteredVel,Time.t(1:end-1));

figure()
subplot(3,1,1)
plot(Time.t(1:end-1),excitationForce,Time.t(1:end-1),mooringForce,Time.t(1:end-1),hydrostaticForce,Time.t(1:end-1),dragForce,Time.t(1:end-1),radForce)
legend('ext','moor','hydr','drag','rad')
ylabel('Sim')
% xlabel('time [s]')
title('Force contributions')

subplot(3,1,2)
plot(Time.t(1:end-1),excitationForceNoisy,Time.t(1:end-1),mooringForceNoisy,Time.t(1:end-1),hydrostaticForceNoisy,Time.t(1:end-1),dragForceNoisy,Time.t(1:end-1),radForceNoisy)
% legend('ext','moor','hydr','drag','rad')
ylabel('Noisy')
% xlabel('time [s]')
% title('Simulated force contributions')

subplot(3,1,3)
plot(Time.t(1:end-1),excitationForceFilter,Time.t(1:end-1),mooringForceFilter,Time.t(1:end-1),hydrostaticForceFilter,Time.t(1:end-1),dragForceFilter,Time.t(1:end-1),radForceFilter)
% legend('ext','moor','hydr','drag','rad')
ylabel('Filtered')
xlabel('time [s]')
% title('Simulated force contributions')


%% Compiling various simulations for plotting in the Paper
% load('EKFSS.mat')






