function [theta_est, epsilon_est] = pf_chem(y, time, sys_specs, a, b, M)

% Length of data
T = length(y);

% Variances
[var, var_A, tp_idx, cut_off, eps_sat, cov_sat] = sys_specs{:};

% Initialize particles
theta_particles = betarnd(1,1,1,M)/2;
epsilon_particles = exprnd(0.1, 1, M);
theta_est(1) = mean(theta_particles);
epsilon_est(1) = mean(epsilon_particles);
ln_w_eps = log(ones(1,M)/M);

for t = 2:T

    % Which region are we in
    r = region(theta_est(t-1), time(t), tp_idx, cut_off);
    mean_eps = {0.5, eps_sat, epsilon_est(t-1), 0.5};

    % Compute proposal parameters
    [alpha, beta] = beta_parameters(theta_particles, a(r), b(r), var, r, cov_sat);

    % Propose particles
    theta_particles = betarnd(alpha, beta)/2;
    epsilon_particles = exprnd(mean_eps{r}, 1,M);

    % Compute epsilon weights
    [w_cov, ln_w_eps, theta_est(t), epsilon_est(t), epsilon_particles] = compute_weights(y(t), epsilon_particles, theta_particles, var_A, M, ln_w_eps);

    % Resample
    idx_cov = datasample(1:M, 1, 'Weights', w_cov);
    theta_sample = theta_particles(idx_cov);

    % Estimate
    %theta_est(t) = mean(theta_particles);
    
end
