function [ force ] = excitationForce_spec(time)
% Excitation force that is a sum of sine waves multiplied by magnitude

global SpectrumData

FX = SpectrumData.Fe;
W = SpectrumData.W;
Mag = SpectrumData.magnitude;
Phi = SpectrumData.phase;
dW = SpectrumData.dW;

for k = 1:length(W)
    IntegralArg(k) = Mag(k)*FX(k)*exp(1i*(W(k)*time+Phi(k)));
end

force = real(dW*trapz(IntegralArg));



end

