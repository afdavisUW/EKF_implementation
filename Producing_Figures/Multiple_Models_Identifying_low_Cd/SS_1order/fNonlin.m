function [ Zdot ] = fNonlin( tau, Z )
% This function is the nonlinear dynamics of the Morison equation
global Func 
% U1 - displacement
% U2 - velocity
% U3 - Acceleration
x1 = Z(1); 
x2 = Z(2);
x3 = Z(3);
x4 = Z(4);


[U1,U2,U3] = particleDynamics_spectrum(tau);


Zdot = Func.f([x1,x2,x3,x4],U1,U2,U3,tau);

end
