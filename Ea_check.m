close all
clear all
clc

% RUN Ea
load temps_info.mat

for n = 1 : length(temps_strings)
    str = join([temps_strings{n}, 'k.mat']);
    load(str)
    k1_adsorb(n) = kOX_est(1);
    k2_adsorb(n) = kOX_est(2);

    k1_desorb(n) = kXO_est(1);
    k2_desorb(n) = kXO_est(2);

end


% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Which data point to exclude
idx = setdiff(1:N, [6]);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;


[Ea4_des, A, Ea_SE, A_SE, ln_k, Rsq_Ea] = get_Ea(k2_desorb(idx), T(idx), R);
[Ea3_des, A, Ea_SE, A_SE, ln_k, Rsq_Ea] = get_Ea(k1_desorb(idx), T(idx), R);

[Ea2_ads, A, Ea_SE, A_SE, ln_k, Rsq_Ea] = get_Ea(k2_adsorb(idx), T(idx), R);
[Ea1_ads, A, Ea_SE, A_SE, ln_k, Rsq_Ea] = get_Ea(k1_adsorb(idx), T(idx), R);


% plot(k2_desorb)

