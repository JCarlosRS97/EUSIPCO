clear, close all
addpath('Functions/')
% file_k1000 = 'Results/Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt5lr5Nrisx100K100000.mat';
% file_k1 = 'Results/Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt5lr5Nrisx100K1.mat';
file_RIS = 'Results/Sweep_dris10/DoF_RIS_MIMO_Nty8Ntz8Nry8Nrz8D10dris8.75lt2lr2Nrisx50K100000.mat';
file_direct = 'Results/Nty8_Ntz8_Nry8_Nrz8_lt2_lr2_Nrisy50_final/Conf_direct_MIMO_Nty8Ntz8Nry8Nrz8D10lt2lr2Nrisx50K100000.mat';
file_direct_RIS = 'Results/Sweep_dris10/DoF_direct_RIS_MIMO_Nty8Ntz8Nry8Nrz8D10dris8.75lt2lr2Nrisx50K100000.mat';

% addpath(dir_files)


% load("Sweep_k/Conf_RIS_MIMO_Nty16Ntz16Nry16Nrz16D5lt5lr5Nrisx50K1000corre.mat")
% plot(10*log10(K_vec), Indir_output_hold), hold on
% plot(10*log10(K_vec), Indir_output_hold_K1000Solution)
% plot(10*log10(K_vec), Indir_output_hold_anom)
% plot(10*log10(K_vec), Indir_output_hold_focusing),hold off
% legend('Optimal', 'LOS', 'Anomalous reflection', 'Focusing function')

figure(1)

plot(-1,-1, 'k-.'), hold on
plot(-1,-1, 'b-o')
plot(-1,-1, 'r-x')
plot(-1,-1, 'g-')
plot(-1,-1, 'b--o')
plot(-1,-1, 'r--x')
plot(-1,-1, 'g--')



%% Reflected channel
load(file_RIS, 'd_ris_vec', 'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal',...
    'Indir_output_K100000_hold', 'Indir_output_K100000_hold_focusing', 'Indir_output_K100000_hold_anom')

figure(1)
plot(d_ris_vec, N_modes_opt, 'b--'),hold on
plot(d_ris_vec(1:2:end), N_modes_opt(1:2:end), 'bo')
% plot(d_ris_vec, N_modes_RIS_LoS, 'r--')
plot(d_ris_vec, N_modes_focus, 'r--')
plot(d_ris_vec(1:2:end), N_modes_focus(1:2:end), 'rx')
plot(d_ris_vec, N_modes_abnormal, 'g--')

figure(2)
plot(d_ris_vec, Indir_output_K100000_hold, 'b--'),hold on
plot(d_ris_vec(1:2:end), Indir_output_K100000_hold(1:2:end), 'bo')
% plot(d_ris_vec, N_modes_RIS_LoS, 'r--')
plot(d_ris_vec, Indir_output_K100000_hold_focusing, 'r--')
plot(d_ris_vec(1:2:end), Indir_output_K100000_hold_focusing(1:2:end), 'rx')
plot(d_ris_vec, Indir_output_K100000_hold_anom, 'g--')

%% Direct channel
load(file_direct,  'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal',...
    'Indir_output_K100000_hold', 'Indir_output_K100000_hold_focusing', 'Indir_output_K100000_hold_anom')
figure(1)
plot(d_ris_vec, N_modes_opt(end)*ones(size(d_ris_vec)), 'k-.')
figure(1)
plot(d_ris_vec, Indir_output_K100000_hold(end)*ones(size(d_ris_vec)), 'k-.')

%% Direct and reflected channel
load(file_direct_RIS, 'd_ris_vec', 'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal',...
    'Indir_output_K100000_hold', 'Indir_output_K100000_hold_focusing', 'Indir_output_K100000_hold_anom')

figure(1)
plot(d_ris_vec, N_modes_opt, 'b-'),hold on
plot(d_ris_vec(1:2:end), N_modes_opt(1:2:end), 'bo')
% plot(d_ris_vec, N_modes_RIS_LoS, 'r--')
plot(d_ris_vec, N_modes_focus, 'r-')
plot(d_ris_vec(1:2:end), N_modes_focus(1:2:end), 'rx')
plot(d_ris_vec, N_modes_abnormal, 'g-')

figure(2)
plot(d_ris_vec, Indir_output_K100000_hold, 'b-'),hold on
plot(d_ris_vec(1:2:end), Indir_output_K100000_hold(1:2:end), 'bo')
% plot(d_ris_vec, N_modes_RIS_LoS, 'r-')
plot(d_ris_vec, Indir_output_K100000_hold_focusing, 'r-')
plot(d_ris_vec(1:2:end), Indir_output_K100000_hold_focusing(1:2:end), 'rx')
plot(d_ris_vec, Indir_output_K100000_hold_anom, 'g-')



%% Tune
figure(1)
legend('Direct channel', 'Direct + reflect channel. Scheme 1.', 'Direct + reflect channel. Scheme 3.',...
    'Direct + reflect channel. Scheme 4.',...
    'Reflect channel. Scheme 1.', 'Reflect channel. Scheme 3.',...
    'Reflect channel. Scheme 4.',...
    'Interpreter', 'latex', 'FontSize', 10)
xlabel('$d_\mathrm{ris}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 12)
axis([1.25,8.75,0.4,4]), grid on