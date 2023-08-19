function [alpha, beta] = beta_parameters(theta_particles, a, b, var, r, cov_sat)

% We are sample 2X ~ Beta(alpha, beta)
if (r == 2)
    mu = 2*normrnd(cov_sat, 0.0001, 1, length(theta_particles));
else  
    mu = 2*(a*theta_particles + b); 
end

if (sum(mu < 0) > 0)
    disp('stop')
end

% Correction test
check = mu - mu.^2;
variances = min( [var*ones(length(theta_particles),1), 0.9*check'.*ones(length(theta_particles),1)]');


alpha = mu.*( mu.*(1 - mu)./variances - 1);
beta = alpha.*(1 - mu)./mu;

if sum( alpha < 0) > 0
    disp('stop')
end

if (sum(isnan(alpha)) > 0 )
    disp('stop')
end

if sum( beta < 0) > 0
    disp('stop')
end

end