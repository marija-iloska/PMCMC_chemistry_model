function [alpha, beta] = beta_parameters(theta_particles, a, b, var, r, cov_sat)

% We are sample 2X ~ Beta(alpha, beta)
if (r == 2)
    mu = 2*normrnd(cov_sat, 0.0001, 1, length(theta_particles));
else
    mu = 2*(a*theta_particles + b);
end

alpha = mu.*( mu.*(1 - mu)./var - 1);
beta = alpha.*(1 - mu)./mu;

end