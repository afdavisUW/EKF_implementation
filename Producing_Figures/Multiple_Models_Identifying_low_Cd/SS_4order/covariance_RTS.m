function Pdot = covariance_RTS( tau, Pvect )
%Function to integrate the covariance backward in time for the RTS filter
% This is equation 5.110 in Cradassis 2nd edition

global Param Time

P = reshape(Pvect, Param.numStates, Param.numStates)';

% Compute F and K
[j_max, k_max] = size(Param.F(:,:,1));
F_current = zeros(j_max*k_max,1);
K_current = zeros(j_max*k_max,1);

count = 1;
for j = 1:j_max
    for k = 1:k_max
        F_temp = squeeze(Param.F(j,k,:));
        K_temp = squeeze(Param.K(j,k,:));
        
        F_current(count) = interp1(Time.hat,F_temp,tau);
        K_current(count) = interp1(Time.hat,K_temp,tau);
        count = count + 1;
    end
end

F = reshape(F_current, Param.numStates, Param.numStates)';
K = reshape(K_current, Param.numStates, Param.numStates)';

% Compute Pdot
PdotMatrix = -(-(F+K)*P - P*(F+K)' + Param.G*Param.Q*Param.G');

Pdot = PdotMatrix(:);

end

