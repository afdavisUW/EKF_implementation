function dxdt = xdot_RTS( tau, x, f )
%Function to integrate the dynamics backwards in time given by equation
% 5.107 in Cradassis 2nd edition This function is used for the RTS filter
global Param Time

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

% Compute U
% [U1,U2,U3] = WaterParticleDyn_vertical(tau, Param.buoyDepth);

% Compute x_f
x_resampled = resample(Param.x.f,tau);
x_f = x_resampled.data;

% Compute particle velocities
horDisp = 0;
[V1,V2,V3,U1,U2,U3] = particleDynamics(horDisp,Time.sim(j));

% Compute xdot
% dxdt = -(F + K)*(x - x_f) - f(x_f);
dxdt = -(-(F + K)*(x - x_f) - f(x_f', U1, U2, U3));

end

