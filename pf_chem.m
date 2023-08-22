function [theta_sample, epsilon_sample] = pf_chem(y, sys_specs, bounds, a, b, M, dat_choice)

% Length of data
T = length(y);

% Variances
[var_A, eps_sat, cov_sat, eps_exp] = sys_specs{:};
[~, cut_off, theta_max, theta_min] = bounds{:};

if (dat_choice == 2)
    rr = [1, 2];
    theta_mean = cov_sat*ones(1,M);
else
    rr = [2, 1];
    theta_mean = 0.1*ones(1,M);
end
r=1;

% Initialize particles
theta_particles = pertrnd(theta_min, theta_mean, theta_max);
epsilon_particles = exprnd(eps_exp, 1, M);
theta_est(1) = mean(theta_particles);
epsilon_est(1) = mean(epsilon_particles);
ln_w_eps = log(ones(1,M)/M);


for t = 2:T

    % Which region are we in
    %r = region(theta_est(t-1), time(t), time(tp_idx), cut_off);
    %mean_eps = {0.5, eps_sat, eps_sat, 0.5};
    mean_eps = {eps_sat, eps_exp};

    % Propose particles
    theta_mean = {cov_sat*ones(1,M), a(r)*theta_particles + b(r)};
    theta_particles = pertrnd(theta_min, theta_mean{rr(r)}, theta_max);
    epsilon_particles = exprnd(mean_eps{rr(r)}, 1,M);

    % Compute epsilon weights
    [w_cov, w_eps, theta_est(t), epsilon_particles, theta_particles] = compute_weights(y(t), epsilon_particles, theta_particles, var_A, M);

    % Store samples
    theta_store(t, :) = theta_particles;
    epsilon_store(t,:) = epsilon_particles;  

    % Estimate
    %theta_est(t) = mean(theta_particles);
    epsilon_est(t) = mean(epsilon_particles.*w_eps);
    idx_eps = datasample(1:M, M, 'Weights', w_eps);
    epsilon_particles = epsilon_particles(idx_eps);

    if dat_choice == 2
        if (theta_est(t) < cut_off)
            r = 2;
        end
    else
        if (theta_est(t) > cut_off)
            r = 2;
        end
    end

    
end

% Take one sample of entire Time horizon
idx = datasample(1:M, 1, 'Weights', w_cov);
theta_sample = theta_store(:, idx);

idx_eps = datasample(1:M, 1, 'Weights', w_eps);
epsilon_sample = epsilon_store(:, idx_eps);

