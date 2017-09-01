function [ Zdot ] = fNonlinIrreg( tau, Z )
% This function is the nonlinear dynamics of the Morison equation
% global Func 
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

[U1,U2,U3] = particleDynamics_spectrum(tau);


Zdot = f_irregular([x1,x2,x3,x4,x5,x6,x7],U1,U2,U3,tau);

end
