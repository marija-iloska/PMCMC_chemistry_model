clear all
close all
clc

% Particle Metropolis within Gibbs
load area_ref490.mat
load epsilons_mat.mat
load temps_info.mat

% System specifications
tp_idx = 45;
cut_off = 0.33;
idx_T = 1;

% Data
time = time_mat_area{idx_T};
T = length(time);
y = area{idx_T};

% Some priors
eps_sat = mean(y(tp_idx - 20 : tp_idx))/cov_sat(idx_T);
eps_exp = [0.5, 0.5, 0.6, 0.6, 0.5, 0.5];

% Number of particles
M = 100;

% Noise
var_A = (std(y(T-10:T)))^2;

% Divide data
y1 = y(1:tp_idx);
y2 = y(tp_idx+1 : T);
time1 = time(1:tp_idx);
time2 = time(tp_idx+1:T);
T1 = length(y1);
T2 = length(y2);


% Bounds for parameter a1 and a2
a_low = 0.001*ones(1,2);
a_high = 0.99*ones(1,2);

% Sample param a1 and a2
a1 = unifrnd(a_low, a_high);

% Bounds for state theta (coverage)
theta_max = 0.5;
theta_min = 0;

% System specifications
sys_specs = {var_A, eps_sat, cov_sat(idx_T), eps_exp(idx_T)};
bounds = {tp_idx, cut_off, theta_max, theta_min};


% Compute bounds for param b
b_low = -a1.*theta_min;
b_high = 0.5 - a1.*theta_max;
mu = (b_high + b_low)/2;

% Sample b
b1 = pertrnd(b_low, mu, b_high);

% Metropolis settings
I = 1;

% Run GIBBS
J = 500;
J0 = round(J/2);

for j = 1:J

    % Sample coverage and absorptivity from PF
    [theta_sample, epsilon_sample] = pf_chem(y1, sys_specs, bounds, a1, b1, M, 1);
    theta1_chain(j,:) = theta_sample;
    epsilon1_chain(j,:) = epsilon_sample;
   
  
    % Sample model parameters a and b for all regions
    tp_AB = find(theta_sample < cut_off);
    tp_AB = tp_AB(end);
    regions = {1 : tp_AB, tp_AB +1 : tp_idx};
    for r = 1:2   
        x = MH_chem(theta_sample, I, r, regions, a1(r), b1(r), a_low(r), a_high(r), bounds, cov_sat(idx_T), 1);
        a1chain(j,r) = x(1);
        b1chain(j,r) = x(2);
        a1(r) = x(1);
        b1(r) = x(2);
    end

    % Posterior for measurement variance
    beta_A = 1./var_A + 0.5*sum( (y1' - theta_sample.*epsilon_sample).^2 );
    var_a1(j) = 1./gamrnd(1, beta_A);
    sys_specs = {var_a1(j), eps_sat, cov_sat(idx_T), eps_exp(idx_T)};
    
    
end

% GET ESTIMATES
a1_est = mean(a1chain(J0:J,:), 1);
b1_est = mean(b1chain(J0:J,:), 1);

% BOUNDARIES FOR NEW PARAMETERS
a_low = a1_est([2, 1]);
a2 = unifrnd(a_low, a_high);

% Compute bounds for param b
b_low = -a2.*theta_min;
b_high = 0.5 - a2.*theta_max;
mu = (b_high + b_low)/2;

% Sample b
b2 = pertrnd(b_low, mu, b_high);


for j = 1:J

    % Sample coverage and absorptivity from PF
    [theta_sample, epsilon_sample] = pf_chem(y2, sys_specs, bounds, a2, b2, M, 2);
    theta2_chain(j,:) = theta_sample;
    epsilon2_chain(j,:) = epsilon_sample;
   
  
    if (j < 2)
    % Sample model parameters a and b for all regions
        tp_AB = find(theta_sample < cut_off);
        tp_AB = tp_AB(1);
        regions = {1 : tp_AB, tp_AB +1 : T2};
    end
    for r = 1:2
    
        x = MH_chem(theta_sample, I, r, regions, a2(r), b2(r), a_low(r), a_high(r), bounds, cov_sat(idx_T), 2);
        a2chain(j,r) = x(1);
        b2chain(j,r) = x(2);
        a2(r) = x(1);
        b2(r) = x(2);
    end

    % Posterior for measurement variance
    beta_A = 1./var_A + 0.5*sum( (y2' - theta_sample.*epsilon_sample).^2 );
    var_a2(j) = 1./gamrnd(1, beta_A);
    sys_specs = {var_a2(j), eps_sat, cov_sat(idx_T), eps_exp(idx_T)};

    % k constants samples
    [c2(j,:), kXO(j,:)] = k_des(a2, b2);
    [c1(j,:), kOX(j,:)] = k_ads(a1_est, b1_est, kXO(j,:));
    
    
end




% % GET ESTIMATES
a2_est = mean(a2chain(J0:J,:), 1);
b2_est = mean(b2chain(J0:J,:), 1);
kXO_est = mean(kXO(J0:J,:), 1);


[c1_est, kOX_est] = k_ads(a1_est, b1_est, kXO_est);


% Get estimates
theta1_est = mean(theta1_chain(J0:J, :),1);
epsilon1_est = mean(epsilon1_chain(J0:J, :),1);

theta2_est = mean(theta2_chain(J0:J, :),1);
epsilon2_est = mean(epsilon2_chain(J0:J, :),1);

theta_est = [theta1_est theta1_est(end), theta2_est(2:end)];
epsilon_est = [epsilon1_est epsilon1_est(end), epsilon2_est(2:end)];


figure;
plot(time1, theta1_chain(1,:), 'Color', 'b', 'LineWidth',2)
hold on
plot(time2(2:end), theta2_chain(1,2:end), 'Color', 'k', 'LineWidth',2)
hold on
plot(time, theta_est, 'Linewidth',2)
title('Coverage', 'FontSize', 15)

figure;
plot(time, epsilon_est.*theta_est)
hold on
plot(time, y)
title('Data vs Estimate', 'FontSize', 15)

figure;
plot(time, epsilon_est)
hold on
plot(time, y)
title('Epsilon', 'FontSize', 15)

figure;
plot(a2chain(:,1), 'linewidth', 1)
hold on
plot(b2chain(:,1), 'linewidth', 1)
title('Chain R3 ', 'FontSize', 15)

figure;
plot(a1chain(:,1), 'linewidth', 1)
hold on
plot(b1chain(:,1), 'linewidth', 1)
title('Chain R1', 'FontSize', 15)



figure;
plot(a2chain(:,2), 'linewidth', 1)
hold on
plot(b2chain(:,2), 'linewidth', 1)
title('Chain R4', 'FontSize', 15)

figure;
plot(a1chain(:,2), 'linewidth', 1)
hold on
plot(b1chain(:,2), 'linewidth', 1)
title('Chain R2', 'FontSize', 15)

figure;
plot(var_a1)
hold on
plot(var_a2)
title('Variance of Measurement Noise', 'FontSize', 15)


figure;
hist(kXO(:,2))
title('R4 Desorption', 'FontSize', 15)

% figure;
% plot(kXO(:,2))


figure;
hist(kXO(:,1))
title('R3 Desorption', 'FontSize', 15)

figure;
hist(kOX(J0:J,1))
title('R1 Adsorption', 'FontSize', 15)


figure;
hist(kOX(J0:J,2))
title('R2 Adsorption', 'FontSize', 15)


% filename = join([temps_strings{idx_T}, 'k.mat']);
% save(filename, 'kOX_est', "kXO_est")

