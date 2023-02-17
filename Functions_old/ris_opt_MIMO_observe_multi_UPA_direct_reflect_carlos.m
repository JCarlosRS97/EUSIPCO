function [ARoutput,AR_K1000_output,AR_anom_output, AR_focusing_output, phase_output, ...
    SingularInfo_OptPhase, SingularInfo_LoSOpt, SingularInfo_abnormal, SingularInfo_focusing, ...
    ObAtRIS_OptPhase, ObAtRIS_LoSOpt,ObAtRIS_abnormal, ObAtRIS_focusing ] ...
    = ris_opt_MIMO_observe_multi_UPA_direct_reflect_carlos(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,SNR,Kdir,K1,K2,f,D,dist_ris,lr,lt,Pt, no_iter,no_mat)
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
    end
end

% Distances
d1 = sqrt(dist_ris^2+lt^2);
d2 = sqrt((D-dist_ris)^2+lr^2);
d_dir = sqrt(D^2+(lt-lr)^2);
%% Scheme 3: RIS optimization
disp('RIS optimization. Scheme 3.')
sin_inc = dist_ris/sqrt( dist_ris^2+lt^2);
sin_ref = -(D-dist_ris)/sqrt( (D-dist_ris)^2+lr^2);

abnormal_phase = -k*(d_dir-d1-d2+(sin_inc+sin_ref)*RISPosition(1,:));

%% Channel generation
disp('Channel generation')
[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_UPA_rec_ob(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,no_mat,Kdir,K1,K2,f,dist_ris,alpha_dir);
%% Scheme 2: RIS optimization
disp('RIS optimization. Scheme 2.')
[Hdirt_K1000,H1t_K1000,H2t_K1000] = chan_mat_RIS_new_model_UPA_rec_ob(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,1,100000,100000,100000, f,dist_ris,alpha_dir);
Hdir_K1000 = Hdirt_K1000{1};
H1_K1000 = H1t_K1000{1};
H2_K1000 = H2t_K1000{1};

Hdir1 = Hdir_K1000;
% Hdir1 = zeros(size(Hdir_K1000)); % Cancel the direct

Qinit = eye(Nt)*(Pt/Nt);
omega_init = ones(1,Nris);
c = 1;
[dC6,stepsize1,stepsize2,~,RIS_phase_opt_K1000] = GPM_Stefan_carlos(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1_K1000*sqrt(SNR),H2_K1000,no_iter,Qinit*c^2,omega_init/c,c);
RIS_phase_K1000 =   angle(RIS_phase_opt_K1000);
% dC6(end)
% RIS_phase_K1000 = -angle(diag(H2_K1000)*H1_K1000);
%% Scheme 4: RIS optimization
disp('RIS optimization. Scheme 4.')
offset = 0;
[Hdirt_siso,H1t_siso,H2t_siso] = chan_mat_RIS_new_model_UPA_rec_ob(1,1,1,1,Nris_x,Nris_y,lt+lambda*offset,lr-lambda*offset,D,1,1000000,1000000,1000000,f,dist_ris,3);
Hdir_siso = Hdirt_siso{1}; H1_siso = H1t_siso{1}; H2_siso = H2t_siso{1};
Phase_RIS_focusing = k*d_dir-angle( H1_siso .* transpose( H2_siso ) );

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
    Hdir1 = Hdir;
%     Hdir1 = zeros(size(Hdir)); % Cancel the direct channel
    [dC6,stepsize1,stepsize2,~,RIS_phase_opt] = GPM_Stefan_carlos(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
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
    omega_init = exp(1j*angle(RIS_phase_opt) );
    c=1;
    [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    AchievableRate_hold = AchievableRate_hold+Cout(end);
    [U,S,V] = svd(Q);
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk} = H1 * U(:,kk) * S(kk,kk);
    end
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk+StreamNum} = diag(exp(1j*angle(RIS_phase_opt) ) ) * H1 * U(:,kk) * S(kk,kk);
    end
    [U,S,V] = svd(H1);
    for kk = 1:StreamNum
        ObAtRIS_OptPhase{kk+StreamNum*2} =  H1 * V(:,kk) ;
    end

    %% Rate scheme 2
    disp('Rate computation. Scheme 2.')
    %     omega_init = exp(1j*RIS_phase_K1000);
    %     [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    %     AchievableRate_K1000_solution_hold = AchievableRate_K1000_solution_hold+dC6(end);
    H_comp = Hdir1 + H2 * diag(exp(1j*RIS_phase_K1000 ) )* H1;
    [U,S,V] = svd(H_comp);
    SingularInfo_LoSOpt{i}{1} = (diag(S))' ;
    SingularInfo_LoSOpt{i}{2} = U;
    SingularInfo_LoSOpt{i}{3} = V;
    omega_init = exp(1j*RIS_phase_K1000);
    [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    [U,S,V] = svd(Q);
    for kk = 1:StreamNum
        ObAtRIS_LoSOpt{kk} = H1 * U(:,kk) * S(kk,kk);
    end
    AchievableRate_K1000_solution_hold = AchievableRate_K1000_solution_hold+Cout(end);

    %% Rate scheme 3
    disp('Rate computation. Scheme 3.')
    RIS_phase_anom = abnormal_phase;
    %     omega_init = exp(1j*RIS_phase_anom);
    %     [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    H_comp = Hdir1 + H2 * diag(exp(1j*RIS_phase_anom ) )* H1;
    [U,S,V] = svd(H_comp);
    SingularInfo_abnormal{i}{1} = (diag(S))' ;
    SingularInfo_abnormal{i}{2} = U;
    SingularInfo_abnormal{i}{3} = V;
    omega_init = exp(1j*RIS_phase_anom);
    [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    [U,S,V] = svd(Q);
    for kk = 1:StreamNum
        ObAtRIS_abnormal{kk} = H1 * U(:,kk) * S(kk,kk);
    end

    AchievableRate_anom_hold = AchievableRate_anom_hold + Cout(end) ;

    %% Rate scheme 4
    disp('Rate computation. Scheme 4.')
    %     omega_init = exp(1j*Phase_RIS_focusing);
    %     [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    %     AchievableRate_focusing_hold = AchievableRate_focusing_hold+dC6(end);
    H_comp = Hdir1 + H2 * diag(exp(1j*Phase_RIS_focusing ) )* H1;
    [U,S,V] = svd(H_comp);
    SingularInfo_focusing{i}{1} = (diag(S))'  ;
    SingularInfo_focusing{i}{2} = U;
    SingularInfo_focusing{i}{3}  = V;
    omega_init = exp(1j*Phase_RIS_focusing);
    [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    [U,S,V] = svd(Q);
    for kk = 1:StreamNum
        ObAtRIS_focusing{kk} = H1 * U(:,kk) * S(kk,kk);
    end
    AchievableRate_focusing_hold = AchievableRate_focusing_hold+Cout(end);

end
phase_output = [ phase_output; transpose( RIS_phase_opt_K1000)];
phase_output = [ phase_output; ( exp(1j*RIS_phase_anom))];

ARoutput = AchievableRate_hold/no_mat;
AR_K1000_output = AchievableRate_K1000_solution_hold/no_mat;
AR_anom_output = AchievableRate_anom_hold/no_mat;
AR_focusing_output = AchievableRate_focusing_hold/ no_mat;