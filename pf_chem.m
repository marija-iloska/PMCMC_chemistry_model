function [theta_est, epsilon_est] = pf_chem(y, time, sys_specs, a, b, M)

% Length of data
T = length(y);

% Variances
[var, var_A, tp_idx, cut_off] = sys_specs{:};

% Initialize particles
theta_particles = betarnd(1,1,1,M)/2;
epsilon_particles = exprnd(0.1, 1, M);
theta_est(1) = mean(theta_particles);
epsilon_est(1) = mean(epsilon_particles);


for t = 2:T

    % Which region are we in
    r = region(theta_est(t-1), time(t), tp_idx, cut_off);

    % Compute proposal parameters
    [alpha, beta] = beta_parameters(theta_particles, a(r), b(r), var);

    % Propose particles
    theta_particles = betarnd(alpha, beta)/2;
    epsilon_particles = exprnd(epsilon_est(t-1), 1,M);

    % Compute epsilon weights
    [w_cov, theta_est(t), epsilon_est(t), epsilon_particles] = compute_weights(y(t), epsilon_particles, theta_particles, var_A, M);

    % Resample
    %idx_cov = datasample(1:M, M, 'Weights', w_cov);
    %theta_particles = theta_particles(idx_cov);

    % Estimate
    %theta_est(t) = mean(theta_particles);
    
end
