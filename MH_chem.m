function [x] = MH_chem(y, I, var, r, regions, a,b, a_low, a_high, cut_off, cov_sat)

% Restore last values
x_old = [a,b];
x = [a; b];
[alpha_old, beta_old] = beta_parameters(y(regions{r}), x_old(1), x_old(2), var, r, cov_sat);

% Max and Min coverages in each region
theta_max = [cut_off, 0.5, 0.5, cut_off];
theta_min = [0, cut_off, cut_off, 0];

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample for a
        a_proposed = unifrnd(a_low, a_high);

        % Compute parameters for pert pdf for b
        b_low = -a_proposed*theta_min(r);
        b_high = 0.5 - a_proposed.*theta_max(r);
        mu = (b_high - b_low)/2;

        % Propose sample for b
        [b_proposed, alpha_pert, beta_pert] = pertrnd( b_low, mu, b_high);
        
        % Store samples
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters
        [alpha_star, beta_star] = beta_parameters(y(regions{r}), x_proposed(1), x_proposed(2), var, r, cov_sat);
        
        
        % LOG Likelihood ratio
        AR1 = sum( (alpha_star - alpha_old).*log(2*y(regions{r})) + (beta_star - beta_old).*log(1 - 2*y(regions{r})));
        
        % LOG Proposal ratio
        AR2 = (alpha_pert - 1)*log((x_proposed(2) - b_low)./(x_old(2) - b_low));
        AR2 = AR2 + (beta_pert - 1)*log( (b_high - x_proposed(2)) ./ (b_high - x_old(2)) );
        
        % Acceptance ratio
        AR = exp(AR1 + AR2);

        % Accept or Reject
        if (rand < AR)
            x(:,i) = x_proposed;
            alpha_old = alpha_star;
            beta_old = beta_star;
        else
            x(:,i) = x_old;
        end

end

end