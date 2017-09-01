function [ Pdot ] = CovarianceSolver( time, Pvec)
% discuss this later

global Fsim Tsim G Q
n = sqrt(length(Pvec));

P = reshape(Pvec, n, n)';

[j_max, k_max] = size(Fsim(:,:,1));
count = 1;
for j = 1:j_max
    for k = 1:k_max
        FsimTemp = squeeze(Fsim(j,k,:));
        F_current(count) = interp1(Tsim,FsimTemp,time);
        count = count + 1;
    end
end

F_current = reshape(F_current, n, n)';

PdotMatrix = F_current*P + P*F_current' + G*Q*G';

Pdot = PdotMatrix(:);

end

