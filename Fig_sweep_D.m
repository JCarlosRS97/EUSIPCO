clear, close all
addpath('Functions/')
dir_files = 'Results/Nty8_Ntz8_Nry8_Nrz8_lt2_lr2_Nrisy50_final/';
% addpath(dir_files)


% load("Sweep_k\Conf_RIS_MIMO_Nty16Ntz16Nry16Nrz16D5lt5lr5Nrisx50K1000corre.mat")
% plot(10*log10(K_vec), Indir_output_hold), hold on
% plot(10*log10(K_vec), Indir_output_hold_K1000Solution)
% plot(10*log10(K_vec), Indir_output_hold_anom)
% plot(10*log10(K_vec), Indir_output_hold_focusing),hold off
% legend('Optimal', 'LOS', 'Anomalous reflection', 'Focusing function')

figure(10)

plot(-1,-1, 'k-.'), hold on
plot(-1,-1, 'b-o')
plot(-1,-1, 'r-x')
plot(-1,-1, 'g-')
plot(-1,-1, 'b--o')
plot(-1,-1, 'r--x')
plot(-1,-1, 'g--')

% plot(-1,-1, 'b-o')
% plot(-1,-1, 'b--o')
% plot(-1,-1, 'r-x')
% plot(-1,-1, 'r--x')
% plot(-1,-1, 'g-x')
% plot(-1,-1, 'g--x')
%% Reflect channel
files = dir(strcat(dir_files, 'Conf_RIS_MIMO*.mat'));

D_vec = cellfun(@(x)str2double(x), regexp( {files.name}, 'D(\d+)', 'tokens', 'once' ));

[D_vec, ind] = sort(D_vec);
files = files(ind);

for i = 1:length(files)
    load(strcat(dir_files, files(i).name), 'Indir_output',...
        'Indir_output_K1000', 'AR_focusing_output', 'AR_anom_output',...
        'SingularInfo_OptPhase_K100000','SingularInfo_LoSOpt_K100000',...
        'SingularInfo_focusing_K100000','SingularInfo_abnormal_K100000');
    rate_RIS_opt(i) = Indir_output;
    rate_RIS_LOS(i) = Indir_output_K1000;
    rate_RIS_focusing(i) = AR_focusing_output;
    rate_RIS_anom(i) = AR_anom_output;

    % Optimal
    N_modes_RIS_opt(i) = effect_rank(SingularInfo_OptPhase_K100000{1}{1}');

    % LOS
    N_modes_RIS_LoS(i) = effect_rank(SingularInfo_LoSOpt_K100000{1}{1}');

    % Focus
    N_modes_RIS_focus(i) = effect_rank(SingularInfo_focusing_K100000{1}{1}');

    % Abnormal
    N_modes_RIS_abnormal(i) = effect_rank(SingularInfo_abnormal_K100000{1}{1}');

end
figure()

plot(D_vec, rate_RIS_opt),hold on
plot(D_vec, rate_RIS_LOS)
plot(D_vec, rate_RIS_focusing)
plot(D_vec, rate_RIS_anom)

legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection')

figure(10)
plot(D_vec, N_modes_RIS_opt, 'b--o'),hold on
% plot(D_vec, N_modes_RIS_LoS, 'r--')
plot(D_vec, N_modes_RIS_focus, 'r--x')
plot(D_vec, N_modes_RIS_abnormal, 'g--')

% legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection',...
%     'Interpreter', 'latex', 'FontSize', 12)
% xlabel('$D$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 14)
% axis([1,10,0,24])

%% Direct channel

files = dir(strcat(dir_files, 'Conf_direct_MIMO*.mat'));

D_vec = cellfun(@(x)str2double(x), regexp( {files.name}, 'D(\d+)', 'tokens', 'once' ));

[D_vec, ind] = sort(D_vec);
files = files(ind);

for i = 1:length(files)
    load(strcat(dir_files, files(i).name), 'Indir_output',...
        'Indir_output_K1000', 'AR_focusing_output', 'AR_anom_output',...
        'SingularInfo_OptPhase_K100000','SingularInfo_LoSOpt_K100000',...
        'SingularInfo_focusing_K100000','SingularInfo_abnormal_K100000');
    rate_direct_opt(i) = Indir_output;
    rate_direct_LOS(i) = Indir_output_K1000;
    rate_direct_focusing(i) = AR_focusing_output;
    rate_direct_anom(i) = AR_anom_output;

    % Optimal
    N_modes_direct_opt(i) = effect_rank(SingularInfo_OptPhase_K100000{1}{1}');

    % LOS
    N_modes_direct_LoS(i) = effect_rank(SingularInfo_LoSOpt_K100000{1}{1}');

    % Focus
    N_modes_direct_focus(i) = effect_rank(SingularInfo_focusing_K100000{1}{1}');

    % Abnormal
    N_modes_direct_abnormal(i) = effect_rank(SingularInfo_abnormal_K100000{1}{1}');

end
figure()

plot(D_vec, rate_direct_opt),hold on
plot(D_vec, rate_direct_LOS)
plot(D_vec, rate_direct_focusing)
plot(D_vec, rate_direct_anom)

legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection')

figure(10)
plot(D_vec, N_modes_direct_opt, 'k-.')
% plot(D_vec, N_modes_direct_LoS, 'r:')
% plot(D_vec, N_modes_direct_focus, 'g:')
% plot(D_vec, N_modes_direct_abnormal, 'm:')

% legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection',...
%     'Interpreter', 'latex', 'FontSize', 12)
% xlabel('$D$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 14)

%% Reflect + Direct channel
files = dir(strcat(dir_files, 'Conf_direct_RIS_MIMO*.mat'));

D_vec = cellfun(@(x)str2double(x), regexp( {files.name}, 'D(\d+)', 'tokens', 'once' ));

[D_vec, ind] = sort(D_vec);
files = files(ind);

for i = 1:length(files)
    load(strcat(dir_files, files(i).name), 'Indir_output',...
        'Indir_output_K1000', 'AR_focusing_output', 'AR_anom_output',...
        'SingularInfo_OptPhase_K100000','SingularInfo_LoSOpt_K100000',...
        'SingularInfo_focusing_K100000','SingularInfo_abnormal_K100000');
    rate_direct_RIS_opt(i) = Indir_output;
    rate_direct_RIS_LOS(i) = Indir_output_K1000;
    rate_direct_RIS_focusing(i) = AR_focusing_output;
    rate_direct_RIS_anom(i) = AR_anom_output;

    % Optimal
    N_modes_direct_RIS_opt(i) = effect_rank(SingularInfo_OptPhase_K100000{1}{1}');

    % LOS
    N_modes_direct_RIS_LoS(i) = effect_rank(SingularInfo_LoSOpt_K100000{1}{1}');

    % Focus
    N_modes_direct_RIS_focus(i) = effect_rank(SingularInfo_focusing_K100000{1}{1}');

    % Abnormal
    N_modes_direct_RIS_abnormal(i) = effect_rank(SingularInfo_abnormal_K100000{1}{1}');

end
figure()

plot(D_vec, rate_direct_RIS_opt),hold on
plot(D_vec, rate_direct_RIS_LOS)
plot(D_vec, rate_direct_RIS_focusing)
plot(D_vec, rate_direct_RIS_anom)

legend('Optimal scheme', 'LOS scheme', 'Focusing scheme', 'Anomolous reflection')

figure(10)
plot(D_vec, N_modes_direct_RIS_opt, 'b-o'),hold on
% plot(D_vec, N_modes_direct_RIS_LoS, 'r-')
plot(D_vec, N_modes_direct_RIS_focus, 'r-x')
plot(D_vec, N_modes_direct_RIS_abnormal, 'g-')

legend('Direct channel', 'Direct + reflect channel. Scheme 1.', 'Direct + reflect channel. Scheme 3.',...
    'Direct + reflect channel. Scheme 4.',...
    'Reflect channel. Scheme 1.', 'Reflect channel. Scheme 3.',...
    'Reflect channel. Scheme 4.',...
    'Interpreter', 'latex', 'FontSize', 10)
xlabel('$D$ (m)', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N_{DoF}$', 'Interpreter', 'latex', 'FontSize', 12)
axis([4,24,0,9]), grid on