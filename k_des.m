function [c, kXO] = k_des(a, b)

% System specifications
dt = 0.067;

% Initialize variables
c = zeros(1,2);
kXO=zeros(1,2);

% SAMPLES_____________________________
% Region III
kXO(1) = (1 - a(1))/dt;
c(1) = b(1)/dt;

% Region IV
kXO(2) = (1 - a(2))/dt;
c(2) = b(2)/dt;






end