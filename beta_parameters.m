function [alpha, beta] = beta_parameters(theta_particles, a, b, var, r, cov_sat)

% We are sample 2X ~ Beta(alpha, beta)
if (r == 2)
    mu = 2*normrnd(cov_sat, 0.0001, 1, length(theta_particles));
else
    mu = 2*(a*theta_particles + b);
end

if (mu < 0)
    disp('stop')
end

% Correction test
check = mu - mu.^2;
if (sum(var > check) > 0)
    var = min(check) - eps
end



alpha = mu.*( mu.*(1 - mu)./var - 1);
beta = alpha.*(1 - mu)./mu;

if sum( alpha < 0) > 0
    disp('stop')
end

if sum( beta < 0) > 0
    disp('stop')
end

end