function [] = readIncidentData(name, samplingFreq, mooringBool)
%readIncidentData is a function that is meant to take in the name of a csv
%dataset from the nnmrec open data project from Brett Bosma
%     This script will read the csv and generate a global variable of
%     sensor data, this will also apply the slope and offset for the load
%     cell datasets
global Data

% displaying correct formats
% mooringBool = [0]; % 0-> taught & 1-> catenary
% samplingFreq.motion = 200;
% samplingFreq.wave = 128;
% name.sensorFile = 'Run_113.csv';
% name.waveFile = 'Run_113.txt';
% name.motionFile = 'Run_113_6D.tsv'; %NOTES: this file needs to be opened and the
% first lines must be deleted. The column at Q needs to be resituated so
% that the same numbers of deliminiters is used on each line CHECK ALSO
% check the number of frames and frequency


%% load cell data from a .csv file
% load load-Cell data needs conversion with sensorCal.mat
load sensorCal.mat % loads the structure of slopes and intercepts for the load cell data
loadCellData = csvread(['dataFiles/',name.sensorFile],1,0); % starts reading at second row, first column
% Data.loadCell.time = loadCellData(:,1); % time data cannot be easily read
% in because the saving format only kept 5 digits of precision, so the time
% string must be re-generated using 128 Hz sampling frequency
Data.loadCell.time = [0: 1/samplingFreq.loadCell : length(loadCellData(:,1))/samplingFreq.loadCell-1/samplingFreq.loadCell]';
Data.loadCell.FPS.port = sensorCal.FPS.port.slope.*loadCellData(:,2)+sensorCal.FPS.port.offset;
Data.loadCell.FPS.bow = sensorCal.FPS.bow.slope.*loadCellData(:,3)+sensorCal.FPS.bow.offset;
Data.loadCell.FPS.stbd = sensorCal.FPS.stbd.slope.*loadCellData(:,4)+sensorCal.FPS.stbd.offset;
Data.loadCell.OWC.bow = sensorCal.OWC.bow.slope.*loadCellData(:,5)+sensorCal.OWC.bow.offset;
Data.loadCell.OWC.port = sensorCal.OWC.port.slope.*loadCellData(:,6)+sensorCal.OWC.port.offset;
Data.loadCell.OWC.stbd = sensorCal.OWC.stbd.slope.*loadCellData(:,7)+sensorCal.OWC.stbd.offset;
Data.loadCell.OWC.umb = sensorCal.FPS.OWCUmb.slope.*loadCellData(:,8)+sensorCal.FPS.OWCUmb.offset;
Data.loadCell.FPS.umb = sensorCal.FPS.Umb.slope.*loadCellData(:,9)+sensorCal.FPS.Umb.offset;
Data.loadCell.pressure = loadCellData(:,10);

% Process net force values acting on FPS
umbAngle = 40; % degrees (this is a best guess)
if mooringBool == 0
    vertAngle = 30; % degrees
else
    vertAngle = 0;
end
horAngle = 90-vertAngle;
angle = 30; % the angle away from the bisecting line between bow and stbd

Data.loadCell.fHeave = cosd(vertAngle).*(Data.loadCell.FPS.bow+Data.loadCell.FPS.port+Data.loadCell.FPS.stbd)+cosd(umbAngle).*Data.loadCell.FPS.umb;
Data.loadCell.fHeaveMean = mean(Data.loadCell.fHeave);
Data.loadCell.fHeave = detrend(Data.loadCell.fHeave,'constant');
Data.loadCell.fHeaveMax = max(abs(Data.loadCell.fHeave));
Data.loadCell.fSurge = cosd(horAngle).*(sind(angle).*(Data.loadCell.FPS.port+Data.loadCell.FPS.stbd)-Data.loadCell.FPS.bow);
Data.loadCell.fSurgeMean = mean(Data.loadCell.fSurge);
Data.loadCell.fSurge = detrend(Data.loadCell.fSurge,'constant');
Data.loadCell.fSurgeMax = max(abs(Data.loadCell.fSurge));

%% motion data frome a .tsv file
motionData = tdfread(['dataFiles/',name.motionFile],'\t'); % residuals and Rot[#]
% values are read into data, but unasigned in the structure currently
noFrames = length(motionData.OWC_X);

% % error check
% if noFrames ~= length(motionData.OWC_X) 
%     display('!ERROR!: number of frames for motion data does not match the number of datapoints for the motion file')
% end

% X -> Surge, Y -> Sway, Z -> Heave
Data.motion.time = [0 : 1/samplingFreq.motion : noFrames/samplingFreq.motion-1/samplingFreq.motion];
Data.motion.OWC.surge = detrend(motionData.OWC_X,'constant')./1000; % /1000 used to convert to meters
Data.motion.OWC.sway = detrend(motionData.Y,'constant')./1000;
Data.motion.OWC.heave = detrend(motionData.Z,'constant')./1000; 
Data.motion.OWC.yaw = detrend(motionData.Yaw,'constant'); % [deg]
Data.motion.OWC.pitch = detrend(motionData.Pitch,'constant'); % [deg]
Data.motion.OWC.roll = detrend(motionData.Roll,'constant'); % [deg]
Data.motion.OWC.residual = motionData.Residual;
Data.motion.FPS.surge = detrend(motionData.FPSPS_X,'constant')./1000;
Data.motion.FPS.surgeMean = mean(motionData.FPSPS_X./1000);
Data.motion.FPS.surgeMax = max(abs(motionData.FPSPS_X./1000));
Data.motion.FPS.sway = detrend(motionData.Y1,'constant')./1000;
Data.motion.FPS.heaveMean = mean(motionData.Z1./1000);
Data.motion.FPS.heaveMax = max(abs(motionData.Z1./1000));
% Data.motion.FPS.heave = detrend(motionData.Z1,'constant')./1000;
Data.motion.FPS.heave = detrend(motionData.Z1)./1000;

Data.motion.FPS.yaw = detrend(motionData.Yaw,'constant');
Data.motion.FPS.pitch = detrend(motionData.Pitch,'constant');
Data.motion.FPS.roll = detrend(motionData.Roll,'constant');
Data.motion.FPS.residual = motionData.Residual1;

span = 0.001;
dt = 1/samplingFreq.motion;
tempFPSHeaveVel = diff(Data.motion.FPS.heave)/dt;
Data.motion.FPS.heaveVel = smooth(Data.motion.time(1:length(tempFPSHeaveVel)),tempFPSHeaveVel,span,'rloess');


%% incident wave data from a text file
fileID = fopen(['dataFiles/',name.waveFile]);
formatSpec = '%s';
N = 17; % number of entries read as stings at beginning of file
waveFileText = textscan(fileID,formatSpec,N,'Delimiter','\t');
waveDataFull = textscan(fileID, '%s %f %f');
% waveData(:,1) = waveDataFull{2}; waveData(:,2) = waveDataFull{3};

Data.wave.time = [0: 1/samplingFreq.wave : length(waveDataFull{2})/samplingFreq.wave-1/samplingFreq.wave]';
Data.wave.incidentOWC = detrend(waveDataFull{2},'constant');
% Data.wave.incidentFPS = detrend(waveDataFull{3},'constant');%+0.007 % adding 7mm
Data.wave.incidentFPS = detrend(waveDataFull{3});


% taking numerical derivitives
span = 0.001; % span percentage for smoothing
dt = 1/samplingFreq.wave;
tempVel = diff(Data.wave.incidentFPS)/dt;
Data.wave.incidentFPSVel = smooth(Data.wave.time(1:length(tempVel)),tempVel,span,'rloess');
% Data.wave.incidentFPSVel = diff(Data.wave.incidentFPS)/dt;

tempAcc = diff(Data.wave.incidentFPSVel)/dt;
Data.wave.incidentFPSAcc = smooth(Data.wave.time(1:length(tempAcc)),tempAcc,span,'rloess');
% Data.wave.incidentFPSAcc = diff(Data.wave.incidentFPSVel)/dt;

% NOTE: attentuation factor still needs to be multiplied in a later
% function

% Generating Bulk parameter for mooring stiffness
% Data.heaveMoorStiff = Data.loadCell.fHeaveMax/Data.motion.FPS.heaveMax;
% Data.surgeMoorStiff = Data.loadCell.fSurgeMax/Data.motion.FPS.surgeMax;
% Data.heaveMoorStiff = Data.loadCell.fHeaveMean/Data.motion.FPS.heaveMean;
% Data.surgeMoorStiff = Data.loadCell.fSurgeMean/Data.motion.FPS.surgeMean;
% for k = 1:length(Data.motion.FPS.heave-5) % to account for 1 timespan being longer
%     Data.stiffVector(k) = interp1(Data.loadCell.time,Data.loadCell.fHeave,Data.motion.time(k))/Data.motion.FPS.heave(k);
% end
save('IncidentWaveData','Data')

end

