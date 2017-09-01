function [ Zdot ] = fNonlinLin( tau, Z )
% This function is the nonlinear dynamics of the Morison equation
global FuncLin
% U1 - displacement
% U2 - velocity
% U3 - Acceleration
x1 = Z(1); 
x2 = Z(2);
x3 = Z(3);

[U1,U2,U3] = particleDynamics_spectrum(tau);


Zdot = FuncLin.f([x1,x2,x3],U1,U2,U3);

end
