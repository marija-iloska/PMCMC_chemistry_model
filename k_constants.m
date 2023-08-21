function [c, kOX, kXO, ads1, ads2] = k_constants(a, b)

% System specifications
dt = 0.067;
P = 0.001;
M = 0.5;

% Initialize variables
c = zeros(1,4);
kOX=zeros(1,2);
kXO=zeros(1,4);

% SAMPLES_____________________________
% Region 4
kXO(4) = (1 - a(4))/dt;
c(4) = b(4)/dt;

kXO(1) = kXO(4);

% Region I
kOX(1) = (1 - a(1) - kXO(1)*dt)/(P*dt);
c(1) = b(1)/dt - kOX(1)*P*M;

% Region III
kXO(3) = (1 - a(3))/dt;
c(3) = b(3)/dt;

kXO(2) = kXO(3);

% Region II
kOX(2) = (1 - a(2) - kXO(2)*dt)/(P*dt);
c(2) = b(2)/dt - kOX(2)*P*M;

ads1 = min([1,  1 - kXO(1)*dt]);
ads2 = min([1,  1 - kXO(2)*dt]);




end