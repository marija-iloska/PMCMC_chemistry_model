function [x] = MH_chem(y, I, var, var_prop, r, mean_in, regions, a,b)

x = mean_in*ones(2,I);
x_old = [a,b];
[alpha_old, beta_old] = beta_parameters(y(regions{r}), x_old(1), x_old(2), var);

theta_r = [0.5, 0.5];
lim = [0.2, 0.2];

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample
        a_proposed = unifrnd(0, lim(r));
        b_proposed = unifrnd(- 0.001*a_proposed, 0.5 - theta_r(r)*a_proposed);
        
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters

        [alpha_star, beta_star] = beta_parameters(y(regions{r}), x_proposed(1), x_proposed(2), var);
        
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