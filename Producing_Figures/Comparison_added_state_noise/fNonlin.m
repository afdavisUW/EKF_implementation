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
x5 = Z(5);
x6 = Z(6);
x7 = Z(7);
x8 = Z(8);
x9 = Z(9);
x10 = Z(10);
x11 = Z(11);
x12 = Z(12);
x13 = Z(13);

[U1,U2,U3] = particleDynamics_spectrum(tau);


Zdot = Func.f([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13],U1,U2,U3,tau);

end