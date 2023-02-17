clear, close all
addpath('Functions\')
% file_k1000 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt5lr5Nrisx100K100000.mat';
% file_k1 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt5lr5Nrisx100K1.mat';
file_RIS = 'Results\Sweep_dris\DoF_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt2lr2Nrisx50K100000.mat';
file_direct = 'Results\Direct_8x8\Conf_direct_MIMO_Nty8Ntz8Nry8Nrz8D5lt5lr5Nrisx50K100000.mat';
file_direct_RIS = 'Results\Sweep_dris\DoF_direct_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5dris3.75lt2lr2Nrisx50K100000.mat';

% addpath(dir_files)


% load("Sweep_k\Conf_RIS_MIMO_Nty16Ntz16Nry16Nrz16D5lt5lr5Nrisx50K1000corre.mat")
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
load(file_RIS, 'd_ris_vec', 'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal')


plot(d_ris_vec, N_modes_opt, 'b--o'),hold on
% plot(d_ris_vec, N_modes_RIS_LoS, 'r--')
plot(d_ris_vec, N_modes_focus, 'r--x')
plot(d_ris_vec, N_modes_abnormal, 'g--')

%% Direct channel
load(file_direct,  'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal')
plot(d_ris_vec, N_modes_opt*ones(size(d_ris_vec)), 'k-.')

%% Reflected channel
load(file_direct_RIS, 'd_ris_vec', 'N_modes_opt', 'N_modes_focus', 'N_modes_abnormal')


plot(d_ris_vec, N_modes_opt, 'b-o'),hold on
% plot(d_ris_vec, N_modes_RIS_LoS, 'r--')
plot(d_ris_vec, N_modes_focus, 'r-x')
plot(d_ris_vec, N_modes_abnormal, 'g-'), hold off


%% Tune

legend('Direct channel', 'Direct + reflect channel. Optimal scheme.', 'Direct + reflect channel. Focusing scheme.',...
    'Direct + reflect channel. Anomalous scheme.',...
    'Reflect channel. Optimal scheme.', 'Reflect channel. Focusing scheme.',...
    'Reflect channel. Anomalous scheme.',...
    'Interpreter', 'latex', 'FontSize', 10)
xlabel('$d_{ris}$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 12)
axis([1.25,3.75,1,12]), grid on