clear all
close all
clc

% Particle Metropolis within Gibbs
load area_ref490.mat

% System specifications
tp_idx = 45;
cut_off = 0.33;

% Data
time = time_mat_area{1};
T = length(time);
y = area{1};

% Number of particles
M = 100;

% Intiialize parameters
a = 0.1;
b = 0.1;


% Noise
var_A = 0.005;
var = 0.001;

% System specifications
sys_specs = {var, var_A, tp_idx, cut_off};

tp_AB = [5, 270];


% b = unifrnd(0, 0.1, 1, 2);
% a = unifrnd(-2*b, 1-2*b);

a = unifrnd(0, 0.2, 1,2);
b = unifrnd(-a*0.001, 0.5 - a*0.5);


% Run PF;
J = 1000;
R = 4;
for j = 1: J
    [theta_est, epsilon_est] = pf_chem(y, time, sys_specs, a, b, M);
    
    
    % Metropolis settings
    I = 1;
    %regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx+1 : tp_AB(2), tp_AB(2)+1 : T };
    regions = {1:tp_idx, tp_idx:T };
    
  
    for r = 1:2
        var_prop = 0.001;
        mean_in = 0.1;
    
        x = MH_chem(theta_est, I, var, var_prop, r, mean_in, regions, a(r), b(r));
        a(r) = mean(x(1, end));
        b(r) = mean(x(2, end));

        achain(j,r) = x(1,end);
        bchain(j,r) = x(2, end);

    end

end

figure;
plot(time, theta_est)
hold on
plot(time,epsilon_est)
hold on
plot(time, y)

figure;
plot(achain(:,1))
hold on
plot(achain(:,2))


figure;
plot(bchain(:,1))
hold on
plot(bchain(:,2))

