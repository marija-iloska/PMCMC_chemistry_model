function [w_eps, theta_est, epsilon_est, epsilon_particles, theta_particles] = compute_weights1(y, epsilon_particles, theta_particles, var_A, M)


% % EPS weights
ln_w_eps = -0.5*log(2*pi*var_A) - 0.5*( (y - epsilon_particles.*theta_particles).^2)/var_A;


% Sample one epsilon particle
w_eps = exp(ln_w_eps - max(ln_w_eps));
w_eps = w_eps./sum(w_eps);


idx_eps = datasample(1:M, M, 'Weights', w_eps);
epsilon_particles = epsilon_particles(idx_eps);
theta_particles = theta_particles(idx_eps);
epsilon_est = mean(epsilon_particles);
theta_est = mean(theta_particles);


end