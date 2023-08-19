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
%y(y<0) = 10e-4;

% Some priors
eps_sat = mean(y(tp_idx - 30 : tp_idx))/cov_sat(1);

% Number of particles
M = 200;


% Noise
var_A = (std(y(T-10:T)))^2;
var = [0.0025, 0.0025, 0.0025, 0.0025];



a_low = 0;
a_high = 1;

theta_max = [0.5, 0.5, 0.5, 0.5];
theta_min = [0, 0.005, 0.005, 0];


% System specifications
sys_specs = {var, var_A, eps_sat, cov_sat(1)};
bounds = {tp_idx, cut_off, theta_max, theta_min};

a = unifrnd(a_low, a_high, 1, 4);
b_low = -a.*theta_min;
b_high = 0.5 - a.*theta_max;
mu = (b_high + b_low)/2;

b = pertrnd(b_low, mu, b_high);



% Metropolis settings
I = 1;
R = 4;


% Run GIBBS
J = 100;
J0 = round(J/2);
for j = 1:J
    [theta_sample, epsilon_est] = pf_chem(y, time, sys_specs, bounds, a, b, M);
    theta_chain(j,:) = theta_sample;
   
  
    % Sample for each region
    tp_AB = find(theta_sample > cut_off);
    if (isempty(tp_AB))
        tp_AB = [5, 80];
    else
        tp_AB = [tp_AB(1), tp_AB(end)-1];
    end
    regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx+1 : tp_AB(2), tp_AB(2)+1 : T-1 };
    for r = 1:R
    
        x = MH_chem(theta_chain(j,:), I, var, r, regions, a(r), b(r), a_low, a_high, bounds, cov_sat(1));
        a(r) = mean(x(1, end));
        b(r) = mean(x(2, end));

        achain(j,r) = x(1,end);
        bchain(j,r) = x(2,end);

    end

end

% Get estimates
theta_est = mean(theta_chain(J0:J, :),1);

figure;
plot(time, theta_est)
hold on
plot(time, theta_chain(J0,:), 'Color', 'b', 'LineWidth',2)
hold on
plot(time, theta_chain(J,:), 'Color', 'k', 'LineWidth',2)
title('Coverage', 'FontSize', 15)

figure;
plot(time, epsilon_est)
hold on
plot(time, y)
title('Epsilon', 'FontSize', 15)

figure;
test = mean(theta_chain(end,:),1);
plot(test)
hold on
plot(movmean(test, 4))

figure;
plot(achain(:,1), 'linewidth', 1)
hold on
plot(achain(:,2), 'linewidth', 1)
hold on
plot(achain(:,3), 'linewidth', 1)
hold on
plot(achain(:,4), 'linewidth', 1)
title('Chain a', 'FontSize', 15)



figure;
plot(bchain(:,1), 'linewidth', 1)
hold on
plot(bchain(:,2), 'linewidth', 1)
hold on
plot(bchain(:,3), 'linewidth', 1)
hold on
plot(bchain(:,4), 'linewidth', 1)
title('Chain b', 'FontSize', 15)

