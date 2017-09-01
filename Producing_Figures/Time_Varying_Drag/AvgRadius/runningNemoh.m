clear all
close all

% num_omega = 25 ;
num_omega = 25 ;

num_omega = 75;
% max_period = 10;
max_period =20;
min_omega = 2*pi/max_period; 
% min_period = .2;
min_period = 0.05;
max_omega = 2*pi/min_period;
% rounding omega to the nearest radian/second
min_omega = round(min_omega,1); max_omega = round(max_omega,1); % Warning!!!
% watch these round functions, because the most common error I am getting
% is when the round function takes one of the input frequencies to be zero
omega = linspace(min_omega,max_omega,num_omega); % the last one that worked
% omega = linspace(5.1,15,3);
min_omega
max_omega

% temp omega
% omega = linspace(4,7,10)

num_angle = 1;
min_angle = 0;
max_angle = 0;
angle = linspace(min_angle,max_angle,num_angle);

depth = 0
[A,B,Fe] = Nemoh(omega,angle,depth);

% AddedMassPlate_inf = A(9,9,end)

save('ABFe_frequency.mat','A','B','Fe','omega')


figure()
B_BEM(1,:) = B(3,3,:);
plot(omega,B_BEM(1,:))

figure()
plot(omega,squeeze(A(3,3,:)))

figure()
plot(omega,abs(Fe(:,3)))