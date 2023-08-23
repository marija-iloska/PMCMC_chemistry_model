close all
clear all
clc

% RUN Ea
load temps_info.mat

for n = 1 : length(temps_strings)
    str = join([temps_strings{n}, 'all.mat']);
    load(str)
    k1_adsorb(n) = kOX_est(1);
    k2_adsorb(n) = kOX_est(2);

    k1_desorb(n) = kXO_est(1);
    k2_desorb(n) = kXO_est(2);

    theta{n} = theta_est;

end


% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Which data point to exclude
idx = setdiff(1:N, [4]);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;


[Ea4_des, A4_des, Ea4_SE, A4_SE, ln_k, Rsq_Ea] = get_Ea(k2_desorb(idx), T(idx), R);
[Ea3_des, A3_des, Ea3_SE, A3_SE, ln_k, Rsq_Ea] = get_Ea(k1_desorb(idx), T(idx), R);

[Ea2_ads, A2_ads, Ea2_SE, A2_SE, ln_k, Rsq_Ea] = get_Ea(k2_adsorb(idx), T(idx), R);
[Ea1_ads, A1_ads, Ea1_SE, A1_SE, ln_k, Rsq_Ea] = get_Ea(k1_adsorb(idx), T(idx), R);


for n = 1:length(idx)
    %plot(time_mat_area{idx(n)},theta{idx(n)}, 'linewidth',2)
    %hold on
    plot(theta{idx(n)}, 'linewidth',2)
    hold on
end
legend(temps_strings{idx})
% save('cov_time_stochastic1.mat', 'theta', 'time_mat_area')

% figure
% plot(1./T(idx), k2_desorb(idx), '.')
