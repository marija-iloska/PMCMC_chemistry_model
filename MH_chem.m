function [x] = MH_chem(y, I, r, regions, a,b, a_low, a_high, bounds, cov_sat, dat_choice)

% Restore last values
x_old = [a,b];
x = [a; b];
[~, ~, theta_max, theta_min] = bounds{:};

if (dat_choice == 2)
    rr = [2,1];
else
    rr = [1,2];
end


mu_a = (a_high + a_low)/2;

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample for a
        
        [a_proposed] = pertrnd(a_low, mu_a, a_high);
        ln_q_star_a  = logpertpdf(a_proposed, a_low, mu_a, a_high);
        ln_q_old_a  = logpertpdf(x_old(1), a_low, mu_a, a_high);
        

        % Compute parameters for pert pdf for b
        b_low = -a_proposed*theta_min;
        b_high = 0.5 - a_proposed.*theta_max;
        mu_b = (b_high + b_low)/2;

        % Propose sample for b
        [b_proposed] = pertrnd(b_low, mu_b, b_high);
        ln_q_star_b  = logpertpdf(b_proposed, b_low, x_old(2), b_high);
        ln_q_old_b  = logpertpdf(x_old(2), b_low, b_proposed, b_high);
        
        % Store samples
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters
        if (0)
            mu_proposed = cov_sat;
        else
            mu_proposed = a_proposed*y(regions{r}(1:end-1)) + b_proposed;
        end
        mu_old = x_old(1)*y(regions{r}(1:end-1)) + x_old(2);
       
        
        % LOG Likelihood ratio
        ln_l_star_theta  = logpertpdf(y(regions{r}(2:end)), theta_min, mu_proposed, theta_max);
        ln_l_old_theta  = logpertpdf(y(regions{r}(2:end)), theta_min, mu_old, theta_max);
        AR1 = sum(ln_l_star_theta - ln_l_old_theta);

        % LOG proposal ratio
        AR2 = ln_q_old_b + ln_q_old_a - ln_q_star_b - ln_q_star_a;
        
        
        % Acceptance ratio
        AR = exp(AR1 + AR2);

        if (isreal(AR) == 0)
            AR = 0;
        end

        % Accept or Reject
        if (rand < AR)
            x(:,i) = x_proposed;
            mu_old = mu_proposed;
        else
            x(:,i) = x_old;
        end

end

end