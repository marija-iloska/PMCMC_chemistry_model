function [x] = MH_chem(y, I, var, r, regions, a,b, a_low, a_high, bounds, cov_sat)

% Restore last values
x_old = [a,b];
x = [a; b];
[alpha_old, beta_old] = beta_parameters(y(regions{r}), x_old(1), x_old(2), var(r), r, cov_sat);
[~, ~, theta_max, theta_min] = bounds{:};

% METROPOLIS HASTINGS
for i = 1 : I

        % Propose sample for a
        [a_proposed, alpha_apert, beta_apert] = pertrnd(a_low, x_old(1), a_high);
        %a_proposed = unifrnd(a_low, a_high);

        % Compute parameters for pert pdf for b
        b_low = -a_proposed*theta_min(r);
        b_high = 0.5 - a_proposed.*theta_max(r);
        mu = (b_high + b_low)/2;

        % Propose sample for b
        [b_proposed, alpha_bpert, beta_bpert] = pertrnd(b_low, mu, b_high);
        
        % Store samples
        x_proposed = [a_proposed, b_proposed];

        % Compute beta parameters
        [alpha_star, beta_star] = beta_parameters(y(regions{r}), x_proposed(1), x_proposed(2), var(r), r, cov_sat);
        
        
        % LOG Likelihood ratio
        AR1 = sum( (alpha_star - alpha_old).*log(2*y(regions{r})) + (beta_star - beta_old).*log(1 - 2*y(regions{r})));
        
        % LOG Proposal ratio
        c = x_old(2) - b_low;
        d = b_high - x_old(2);
        if ( c < 0 )
            c = 10e-4;
        end
        if ( d < 0 )
            d = 10e-4;
        end

        AR2 = (alpha_bpert - 1)*log(c./(x_proposed(2) - b_low));
        AR2 = AR2 + (beta_bpert - 1)*log( d./(b_high - x_proposed(2)));

        AR3 = (alpha_apert - 1)*log((x_old(1) - a_low)./(x_proposed(2) - a_low));
        AR3 = AR3 + (beta_apert - 1)*log( (a_high - x_old(1))./(a_high - x_proposed(1)));

% 
%         if (isreal(AR2) == 0)
%             disp('stop')
%         end
        
        % Acceptance ratio
        AR = exp(AR1 + AR2 + AR3);

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