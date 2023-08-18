function [x] = MH_chem(y, I, var, r, regions, a,b, a_low, a_high, cut_off, cov_sat)


x_old = [a,b];
x = [a; b];
[alpha_old, beta_old] = beta_parameters(y(regions{r}), x_old(1), x_old(2), var, r, cov_sat);

theta_r = [cut_off, 0.5, 0.5, cut_off];
lim = [0, cut_off, cut_off, 0];

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample
        a_proposed = unifrnd(a_low, a_high);
        b_proposed = unifrnd(- lim(r)*a_proposed, theta_r(r) - a_proposed*theta_r(r));
        
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters
        [alpha_star, beta_star] = beta_parameters(y(regions{r}), x_proposed(1), x_proposed(2), var, r, cov_sat);
        
        % Acceptance ratio
        AR = sum( (alpha_star - alpha_old).*log(2*y(regions{r})) + (beta_star - beta_old).*log(1 - 2*y(regions{r})));
        AR = exp(AR);

        if (rand < AR)
            x(:,i) = x_proposed;
            alpha_old = alpha_star;
            beta_old = beta_star;
        else
            x(:,i) = x_old;
        end

end

end