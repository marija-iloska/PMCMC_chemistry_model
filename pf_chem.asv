function [theta_sample, epsilon_est] = pf_chem(y, time, sys_specs, bounds, a, b, M)

% Length of data
T = length(y);

% Variances
[var, var_A, eps_sat, cov_sat] = sys_specs{:};
[tp_idx, cut_off, theta_max, theta_min] = bounds{:};

% Initialize particles
theta_particles = betarnd(0.01,1,1,M)/2;
epsilon_particles = exprnd(0.5, 1, M);
theta_est(1) = mean(theta_particles);
epsilon_est(1) = mean(epsilon_particles);
ln_w_eps = log(ones(1,M)/M);

% theta_low = [0, cut_off, cut_off, 0];
% theta_high =[.0, 0.5, 0.5, cut_off];

for t = 2:T

    % Which region are we in
    r = region(theta_est(t-1), time(t), tp_idx, cut_off);
    mean_eps = {0.5, eps_sat, epsilon_est(t-1), 0.5};

    % Compute proposal parameters
    %[alpha, beta] = beta_parameters(theta_particles, a(r), b(r), var(r), r, cov_sat);

    % Propose particles
    theta_mean = a(r)*theta_particles + b(r);
    if (r == 2)
        theta_mean = cov_sat*ones(1,M);
    end
    theta_particles = pertrnd(theta_min(r), theta_mean, theta_max(r));
    %theta_particles = betarnd(alpha, beta)/2;
    epsilon_particles = exprnd(mean_eps{r}, 1,M);

    if ( sum(isnan(theta_particles)) > 0)
        disp('stop')
    end

    % Compute epsilon weights
    [w_cov, ln_w_eps, theta_est(t), epsilon_est(t), epsilon_particles, theta_particles] = compute_weights(y(t), epsilon_particles, theta_particles, var_A, M, ln_w_eps);

    % Resample
    %idx_cov = datasample(1:M, 1, 'Weights', w_cov);
    %theta_particles = theta_particles(idx_cov);
    theta_store(t, :) = theta_particles;

    % Estimate
    %theta_est(t) = mean(theta_particles);
    
end

idx = datasample(1:M, 1, 'Weights', w_cov);
theta_sample = theta_store(:, idx);
