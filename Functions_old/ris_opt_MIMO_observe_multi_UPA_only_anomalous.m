function [ARoutput,AR_K1000_output,AR_anom_output, AR_focusing_output, phase_output, ...
    SingularInfo_abnormal] ...
    = ris_opt_MIMO_observe_multi_UPA_only_anomalous(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,SNR,Kdir,K1,K2,f,D,dist_ris,lr,lt,Pt, no_iter,no_mat)
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
%% Scheme 3: RIS optimization
disp('RIS optimization. Scheme 3.')
sin_inc = dist_ris/sqrt( dist_ris^2+lt^2);
sin_ref = -(D-dist_ris)/sqrt( (D-dist_ris)^2+lr^2);

abnormal_phase = -k*(sin_inc+sin_ref)*RISPosition(1,:);

%% Channel generation
disp('Channel generation')
[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_UPA_rec_ob(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,no_mat,Kdir,K1,K2,f,dist_ris,alpha_dir);

%% Rate computation
AchievableRate_anom_hold = 0;
AchievableRate_K1000_solution_hold = 0;
AchievableRate_focusing_hold = 0;
AchievableRate_hold = 0;

phase_output = [ ];
SingleValue_output = [ ];
    Qinit = eye(Nt)*(Pt/Nt);

for i = 1:no_mat
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
    Hdir1 = zeros(size(Hdir));
    c =1;
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

    AchievableRate_anom_hold = AchievableRate_anom_hold + Cout(end) ;


end

ARoutput = AchievableRate_hold/no_mat;
AR_K1000_output = AchievableRate_K1000_solution_hold/no_mat;
AR_anom_output = AchievableRate_anom_hold/no_mat;
AR_focusing_output = AchievableRate_focusing_hold/ no_mat;