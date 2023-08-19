function [x] = MH_chem(y, I, var, r, regions, a,b, a_low, a_high, bounds, cov_sat)

% Restore last values
x_old = [a,b];
x = [a; b];
%[alpha_old, beta_old] = beta_parameters(y(regions{r}), x_old(1), x_old(2), var(r), r, cov_sat);
[~, ~, theta_max, theta_min] = bounds{:};

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample for a
        [a_proposed] = pertrnd(a_low, x_old(1), 0.99);
        ln_q_star_a  = logpertpdf(a_proposed, a_low, x_old(1), a_high);
        ln_q_old_a  = logpertpdf(x_old(1), a_low, a_proposed, a_high);
        

        % Compute parameters for pert pdf for b
        b_low = -a_proposed*theta_min(r);
        b_high = 0.5 - a_proposed.*theta_max(r);
        mu = (b_high + b_low)/2;

        % Propose sample for b
        [b_proposed] = pertrnd(b_low, mu, b_high);
        ln_q_star_b  = logpertpdf(b_proposed, b_low, x_old(2), b_high);
        ln_q_old_b  = logpertpdf(x_old(2), b_low, b_proposed, b_high);
        
        % Store samples
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters
        %[alpha_star, beta_star] = beta_parameters(y(regions{r}), x_proposed(1), x_proposed(2), var(r), r, cov_sat);
        mu_proposed = a_proposed*y(regions{r}(1:end-1)) + b_proposed;
        mu_old = x_old(1)*y(regions{r}(1:end-1)) + x_old(2);
       
        
        % LOG Likelihood ratio
        %AR1 = sum( (alpha_star - alpha_old).*log(2*y(regions{r})) + (beta_star - beta_old).*log(1 - 2*y(regions{r})));
        ln_l_star_theta  = logpertpdf(y(regions{r}(2:end)), theta_min(r), mu_proposed, theta_max(r));
        ln_l_old_theta  = logpertpdf(y(regions{r}(2:end)), theta_min(r), mu_old, theta_max(r));
        AR1 = sum(ln_l_star_theta - ln_l_old_theta);

        % LOG proposal ratio
        AR2 = ln_q_old_b + ln_q_old_a - ln_q_star_b - ln_q_star_a;
        if (isreal(AR2) == 0)
            AR2 = 0;
        end
        
        % Acceptance ratio
        AR = exp(AR1 + AR2);

        % Accept or Reject
        if (rand < AR)
            x(:,i) = x_proposed;
            mu_old = mu_proposed;
%             alpha_old = alpha_star;
%             beta_old = beta_star;
        else
            x(:,i) = x_old;
        end

end

end