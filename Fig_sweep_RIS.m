clear, close all
addpath('Functions\')
% file_k1000 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt2lr2Nrisx100K100000_final.mat';
% file_k1 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D5lt2lr2Nrisx100K1_final.mat';
file_k1000 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D15lt2lr2Nrisx100K1000000corre.mat';
file_k1 = 'Results\Conf_RIS_MIMO_Nty8Ntz8Nry8Nrz8D15lt2lr2Nrisx100K1corre.mat';
% addpath(dir_files)


% load("Sweep_k\Conf_RIS_MIMO_Nty16Ntz16Nry16Nrz16D5lt5lr5Nrisx50K1000corre.mat")
% plot(10*log10(K_vec), Indir_output_hold), hold on
% plot(10*log10(K_vec), Indir_output_hold_K1000Solution)
% plot(10*log10(K_vec), Indir_output_hold_anom)
% plot(10*log10(K_vec), Indir_output_hold_focusing),hold off
% legend('Optimal', 'LOS', 'Anomalous reflection', 'Focusing function')

figure(1)

plot(-1,-1, 'k-'), hold on
plot(-1,-1, 'k--')
plot(-1,-1, 'b-x')
plot(-1,-1, 'r-o')
plot(-1,-1, 'g-diamond')
plot(-1,-1, 'm-')

% plot(-1,-1, 'b-o')
% plot(-1,-1, 'b--o')
% plot(-1,-1, 'r-x')
% plot(-1,-1, 'r--x')
% plot(-1,-1, 'm-x')
% plot(-1,-1, 'm--x')
%% K10000 channel
    load(file_k1000, 'Nris_y', 'NRIS_oneside','Indir_output_K100000_hold',...
        'Indir_output_K100000_hold_K1000Solution', 'Indir_output_K100000_hold_focusing', 'Indir_output_K100000_hold_anom');
    rate_k1000_opt = Indir_output_K100000_hold;
    rate_k1000_LOS= Indir_output_K100000_hold_K1000Solution;
    rate_k1000_focusing = Indir_output_K100000_hold_focusing;
    rate_k1000_anom = Indir_output_K100000_hold_anom;


plot(NRIS_oneside*Nris_y, rate_k1000_opt,'b--x'),hold on
plot(NRIS_oneside*Nris_y, rate_k1000_LOS,'r--o')
plot(NRIS_oneside*Nris_y, rate_k1000_focusing,'g--diamond')
plot(NRIS_oneside*Nris_y, rate_k1000_anom,'m--')


%% K1 channel

load(file_k1, 'Nris_y', 'NRIS_oneside','Indir_output_K100000_hold',...
        'Indir_output_K100000_hold_K1000Solution', 'Indir_output_K100000_hold_focusing', 'Indir_output_K100000_hold_anom');

    rate_k1_opt = Indir_output_K100000_hold;
    rate_k1_LOS = Indir_output_K100000_hold_K1000Solution;
    rate_k1_focusing = Indir_output_K100000_hold_focusing;
    rate_k1_anom = Indir_output_K100000_hold_anom;


plot(NRIS_oneside*Nris_y, rate_k1_opt,'b-x'),hold on
plot(NRIS_oneside*Nris_y, rate_k1_LOS,'r-o')
plot(NRIS_oneside*Nris_y, rate_k1_focusing,'g-diamond')
plot(NRIS_oneside*Nris_y, rate_k1_anom,'m-')

legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection')

% figure(10)
% plot(NRIS_oneside*Nris_y, N_modes_direct_opt, 'k:')
% plot(NRIS_oneside*Nris_y, N_modes_direct_LoS, 'r:')
% plot(NRIS_oneside*Nris_y, N_modes_direct_focus, 'g:')
% plot(NRIS_oneside*Nris_y, N_modes_direct_abnormal, 'm:')

% legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection',...
%     'Interpreter', 'latex', 'FontSize', 12)
% xlabel('$D$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 14)


%% Tune

legend('K = 1', 'K = 100000','Scheme 1', 'Scheme 2',...
    'Scheme 3','Scheme 4',...
    'Interpreter', 'latex', 'FontSize', 10, 'Location','Northwest')
xlabel('Number of RIS elements', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Capacity (bps)', 'Interpreter', 'latex', 'FontSize', 12)
axis([2,2500, 0,100]),grid on