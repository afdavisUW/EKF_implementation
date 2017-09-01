clear all
% close all
% clc
% help axiMesh
format compact
set(0,'DefaultFigurePosition','factory')

% This code is meant to experiment with creating a mesh

% number = 30;

% R = 0.2706;
Rb = (.2706+0.25)/2;
Rs = (.1894+.175)/2;

% z = [0,-.075,-.075]
% r = [Rb,Rs,0]
z = [0.15,0.15,0,-0.07,-0.075,-0.075];
r = [0,Rb,Rb,Rs+(Rb-Rs)/3,Rs,0];
% z = -linspace(-Rb,Rb,number)
% r = sqrt(Rb^2-z.^2)

[Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,length(z))


% clear all
% % close all
% % clc
% % help axiMesh
% format compact
% set(0,'DefaultFigurePosition','factory')
% 
% % This code is meant to experiment with creating a mesh
% 
% number = 30;
% 
% % z = -linspace(-0.225,0.225,number)
% % r = sqrt(0.225^2-z.^2)
% z = -linspace(-0.2706,.2706,number)
% r = sqrt(.2706^2-z.^2)
% 
% [Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,number)
