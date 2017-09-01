function [U1,U2,U3] = particleDynamics(time)
%This function is meant to take into account the celerity of the wave so
%that based on a given horizontal position you can determine where the
%object is on the wave and therefore the horizontal and vertical dynamics
%of the water particles

% Note positive horizontal displacement is to the right, positive vertical
% displacement is to the left

global Param Data


    % vertical wave dynamics
U1 = interp1(Data.wave.time,Data.wave.incidentFPS,time);
% U2 = Param.attenFactor*interp1(Data.wave.time(1:length(Data.wave.incidentFPSVel)),Data.wave.incidentFPSVel,time);
% U3 = Param.attenFactor*interp1(Data.wave.time(1:length(Data.wave.incidentFPSAcc)),Data.wave.incidentFPSAcc,time);
U2 = interp1(Data.wave.time(1:length(Data.wave.incidentFPSVel)),Data.wave.incidentFPSVel,time);
U3 = interp1(Data.wave.time(1:length(Data.wave.incidentFPSAcc)),Data.wave.incidentFPSAcc,time);



end

