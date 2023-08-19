clear all
close all
clc

% Particle Metropolis within Gibbs
load area_ref490.mat
load expected_coverage.mat

% System specifications
tp_idx = 45;
cut_off = 0.33;

% Data
time = time_mat_area{1};
T = length(time);
y = area{1};
eps_sat = mean(y(tp_idx - 30 : tp_idx))/cov_sat(1);

% Number of particles
M = 500;


% Noise
var_A = 0.01;
var = 0.005;

% System specifications
sys_specs = {var, var_A, tp_idx, cut_off, eps_sat, cov_sat(1)};

tp_AB = [4, 270];


a_low = 0;
a_high = 0.5;

a = unifrnd(a_low, a_high, 1, 4);
b = unifrnd(0, cut_off*(1 - a));
b([2,3]) = unifrnd(-a([2,3])*cut_off, 0.5 - a([2,3])*0.5);


% Run PF;
J = 100;
R = 4;
for j = 1:J
    [theta_est(j,:), epsilon_est] = pf_chem(y, time, sys_specs, a, b, M);
    
    
    % Metropolis settings
    I = 1;
    regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx+1 : tp_AB(2), tp_AB(2)+1 : T-1 };
    
  
    for r = 1:4
    
        x = MH_chem(theta_est, I, var, r, regions, a(r), b(r), a_low, a_high, cut_off, cov_sat(1));
        a(r) = mean(x(1, end));
        b(r) = mean(x(2, end));

        achain(j,r) = x(1,end);
        bchain(j,r) = x(2, end);

    end

end



figure;
plot(time, theta_est)
hold on
plot(time, theta_est(1,:), 'Color', 'b', 'LineWidth',2)
hold on
plot(time, theta_est(J,:), 'Color', 'k', 'LineWidth',2)
plot(time, epsilon_est)
hold on
plot(time, y)

figure;
plot(achain(:,1))
hold on
plot(achain(:,2))
hold on
plot(achain(:,3))
hold on
plot(achain(:,4))
hold on



figure;
plot(bchain(:,1))
hold on
plot(bchain(:,2))
hold on
plot(bchain(:,3))
hold on
plot(bchain(:,4))

