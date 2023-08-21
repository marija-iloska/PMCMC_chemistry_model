function [c, kOX] = k_ads(a, b, kXO)

% System specifications
dt = 0.067;
P = 0.001;
M = 0.5;

% Initialize variables
c = zeros(1,2);
kOX=zeros(1,2);

% SAMPLES_____________________________
% Region I
kOX(1) = (1 - a(1) - kXO(2)*dt)/(P*dt);
c(1) = b(1)/dt - kOX(1)*P*M;

% Region II
kOX(2) = (1 - a(2) - kXO(1)*dt)/(P*dt);
c(2) = b(2)/dt - kOX(2)*P*M;





end