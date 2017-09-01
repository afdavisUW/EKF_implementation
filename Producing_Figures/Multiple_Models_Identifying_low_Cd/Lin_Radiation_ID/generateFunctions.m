function [ func ] = generateFunctions( )
% generate functions populates the working foulder with functions
%   This should be 1 of 2 places that the function needs to be specified.
%   The other location is the nonlinear function that is used in the ode
%   solver.
global Param

%% Define Input and Output Functions
syms symT x1 x2 x3 U1 U2 U3 time
vars = [x1,x2,x3];

% Input and Output function needs to be written here
dynamics = [x2;
    1/Param.totalInertia*(Excite_forceSpectrum(time) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*x1 - Param.linRad*x2 - Param.dragFactor*x3*(x2-U2)*abs(x2-U2));
    0];



output = [x1;x2];

%% Generate .m files for: f, F, h, and H
Jacob_f = jacobian(dynamics,vars);
Jacob_h = jacobian(output,vars);

func.f = matlabFunction(dynamics,'File','f','Vars',{vars,U1,U2,U3,time});
func.h = matlabFunction(output,'File','h','Vars',{vars});

func.F = matlabFunction(Jacob_f,'File','Jacob_f','Vars',{vars,U1,U2,U3});
func.H = matlabFunction(Jacob_h,'File','Jacob_h','Vars',{vars});

end

