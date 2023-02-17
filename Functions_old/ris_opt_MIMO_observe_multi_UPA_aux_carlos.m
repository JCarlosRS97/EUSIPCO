function [ARoutput,AR_K1000_output,AR_anom_output, AR_focusing_output, phase_output, ...
    SingularInfo_OptPhase, SingularInfo_LoSOpt, SingularInfo_abnormal, SingularInfo_focusing, ...
    ObAtRIS_OptPhase, ObAtRIS_LoSOpt,ObAtRIS_abnormal, ObAtRIS_focusing ] ...
    = ris_opt_MIMO_observe_multi_UPA_aux_carlos(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,P_noise,Kdir,K1,K2,f,D,dist_ris,lr,lt,Pt, no_iter,no_mat)
% clear; close all;
% Nt = 8;
% Nr = 4;
% Nris = 15^2;
% SNR = db2pow(120);
% K = 1;
% Nt=1;
% Nr=1;

nonDIR = 0;
% no_mat = 1;
quant_bits = 0;
% no_iter = 500;
alpha_dir = 3;

lambda = 3e8/f;     % wavelength
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber
Nt = Nt_y * Nt_z;
Nr = Nr_y * Nr_z;
Nris = Nris_x*Nris_y;
RISPosition = zeros(3,Nris);
N1 = sqrt(Nris);
for i = 1:Nris_x
    for ii = 1:Nris_y
        RISPosition( 1, (i-1)*Nris_y+ii ) = (i - (Nris_x+1)/2)*dris;
        RISPosition( 2, (i-1)*Nris_y+ii ) = (ii - (Nris_y+1)/2)*dris;
        RISPosition( 2, (i-1)*Nris_y+ii ) = 0;
    end
end

%% Channel generation
disp('Channel generation')
[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_UPA_rec_ob(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,no_mat,Kdir,K1,K2,f,dist_ris,alpha_dir);

%% Focusing scheme 2d: RIS optimization
disp('Focusing scheme 2d.')
offset = 0;
[Hdirt_siso,H1t_siso,H2t_siso] = chan_mat_RIS_new_model_UPA_rec_ob(1,1,1,1,Nris_x,Nris_y,lt+lambda*offset,lr-lambda*offset,D,1,1000000,1000000,1000000,f,dist_ris,3);
Hdir_siso = Hdirt_siso{1}; H1_siso = H1t_siso{1}; H2_siso = H2t_siso{1};
Phase_RIS_focusing = -angle( H1_siso .* transpose( H2_siso ) );

%% Focusing scheme ideal.
center_1 = [-dist_ris, 0, lt];
center_2 = [D-dist_ris, 0, lr];
focusing_f_1 = compute_focusing_function(RISPosition, center_1, k);
focusing_f_2 = compute_focusing_function(RISPosition, center_2, k);

%% Rate computation
AchievableRate_anom_hold = 0;
AchievableRate_K1000_solution_hold = 0;
AchievableRate_focusing_hold = 0;
AchievableRate_hold = 0;

phase_output = [ ];
SingleValue_output = [ ];

for i = 1:no_mat
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};

    %% Scheme 1: RIS optimization
    disp('RIS optimization. Scheme 1.')
    StreamNum = min(Nt,Nr);
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris);
    c = 1;
    Hdir1 = zeros(size(Hdir));
    [dC6,stepsize1,stepsize2,~,RIS_phase_opt] = GPM_Stefan_carlos(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(1/P_noise)/c,H1*sqrt(1/P_noise),H2,no_iter,Qinit*c^2,omega_init/c,c);
    %     AchievableRate_hold = AchievableRate_hold+dC6(end);
    %     AchievableRate_hold = AchievableRate_hold+dC6;
    phase_output = [ phase_output; transpose( RIS_phase_opt)];

    %% Rate scheme 1
    disp('Rate computation. Scheme 1.')
    H_comp = Hdir1 + H2 * diag(exp(1j*angle(RIS_phase_opt) ) )* H1;
    [U,S,V] = svd(H_comp);
    SingularInfo_OptPhase{i}{1} = (diag(S))';
    SingularInfo_OptPhase{i}{2} = U;
    SingularInfo_OptPhase{i}{3} = V;

    % original code
%     omega_init = exp(1j*angle(RIS_phase_opt) );
%     c=1;
%     [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(1/P_noise)/c,H1*sqrt(1/P_noise),H2,no_iter,Qinit*c^2,omega_init/c,c);
%     AchievableRate_hold = AchievableRate_hold+Cout(end);
%     [U,S,V] = svd(Q);

    % Water filling
    D = water_fill(Pt, diag(S).^2/P_noise);
    Q = V*diag(D)*V';
    AchievableRate = real(log(det(eye(Nr)+H_comp*Q*H_comp'/P_noise)))/log(2);
    AchievableRate_hold = AchievableRate_hold + AchievableRate;

    % Modification
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk} = H1 * V(:,kk) * sqrt(D(kk));
    end
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk+StreamNum} = diag(exp(1j*angle(RIS_phase_opt) ) ) * H1 * V(:,kk) * sqrt(D(kk));
    end
    [U,S,V] = svd(H1);
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk+StreamNum*2} =  H1 * V(:,kk) ;
    end


    %% Focusing scheme
    disp('Rate computation. Focusing scheme.')
    %     omega_init = exp(1j*Phase_RIS_focusing);
    %     [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    %     AchievableRate_focusing_hold = AchievableRate_focusing_hold+dC6(end);
    H_comp = Hdir1 + H2 * diag(exp(1j*Phase_RIS_focusing ) )* H1;
    [U,S,V] = svd(H_comp);
    SingularInfo_focusing{i}{1} = (diag(S))'  ;
    SingularInfo_focusing{i}{2} = U;
    SingularInfo_focusing{i}{3}  = V;

    % Original code
%     omega_init = exp(1j*Phase_RIS_focusing);
%     [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(1/P_noise)/c,H1*sqrt(1/P_noise),H2,no_iter,Qinit*c^2,omega_init/c,c);
%     [U,S,V] = svd(Q);

    % Waterfilling
    D = water_fill(Pt, diag(S).^2/P_noise);
    Q = V*diag(D)*V';
    AchievableRate = real(log(det(eye(Nr)+H_comp*Q*H_comp'/P_noise)))/log(2);
    AchievableRate_focusing_hold = AchievableRate_focusing_hold + AchievableRate;

    for kk = 1:StreamNum
        ObAtRIS_focusing{kk} = H1 * V(:,kk) * D(kk);
    end
%     AchievableRate_focusing_hold = AchievableRate_focusing_hold+Cout(end);

end


phase_output = [ phase_output; transpose( RIS_phase_opt_K1000)];
phase_output = [ phase_output; ( exp(1j*RIS_phase_anom))];

ARoutput = AchievableRate_hold/no_mat;
AR_K1000_output = AchievableRate_K1000_solution_hold/no_mat;
AR_anom_output = AchievableRate_anom_hold/no_mat;
AR_focusing_output = AchievableRate_focusing_hold/ no_mat;