function [Pdot] = covariance_EKF( tau, Pvect )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Param Time

P = reshape(Pvect, Param.numStates, Param.numStates)';

[j_max, k_max] = size(Param.Fsim(:,:,1));
% F_current = zeros(j_max*k_max,1);

count = 1;
for j = 1:j_max
    for k = 1:k_max
        FsimTemp = squeeze(Param.Fsim(j,k,:));
        F_current(count) = interp1(Time.sim,FsimTemp,tau);
        count = count + 1;
    end
end

F_current = reshape(F_current, Param.numStates, Param.numStates)';

PdotMatrix = F_current*P + P*F_current' + Param.GQGt;
% PdotMatrix = F_current*P + P*F_current' + Param.G*Param.Q*Param.G';

Pdot = PdotMatrix(:);

end

