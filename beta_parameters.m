function [alpha, beta] = beta_parameters(theta_particles, a, b, var)

% We are sample 2X ~ Beta(alpha, beta)
mu = 2*(a*theta_particles + b);

alpha = mu.*( mu.*(1 - mu)./var - 1);
beta = alpha.*(1 - mu)./mu;

end