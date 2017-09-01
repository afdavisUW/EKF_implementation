function [ Force ] = Excite_forceSpectrum( time  )
%Excitation force funciton
% [ Force ] = Excite_Force(Height, Frequency, time  )

global SpectrumData
FX = SpectrumData.Fe;
W = SpectrumData.W;
Mag = SpectrumData.magnitude;
Phi = SpectrumData.phase;
dW = SpectrumData.dW;

% linear version
% Force = real(Height/2*ExcitationW_vec(1,index)*exp(1i*Frequency*time));

for k = 1:length(W)
%     InsideIntegral(k) = sqrt(2*Mag(k))*FX(k)*exp(1i*(W(k)*time+Phi(k)));
%     The second attempt is actually balancing the units and taking into
%     account the refomation of teh GlobalMagnitude vector as a series of
%     wave amplitudes, rather than energies
%     InsideIntegral(k) = sqrt(2*Mag(k))*FX(k)*exp(1i*(W(k)*time+Phi(k)));
%     InsideIntegral(k) = Mag(k)*FX(k)*exp(1i*(W(k)*time+Phi(k)));
    InsideIntegral(k) = sqrt(2*Mag(k)*dW)*FX(k)*exp(1i*(W(k)*time+Phi(k)));

    
end

% Force = real(dW*trapz(InsideIntegral));
Force = real(trapz(InsideIntegral));


end

