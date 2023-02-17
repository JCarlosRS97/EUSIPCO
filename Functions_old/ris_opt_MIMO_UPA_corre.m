function [RIS_config_OptPhase, RIS_config_LoSOpt, RIS_config_abnormal, RIS_config_focusing, H1,H2  ] ...
    = ris_opt_MIMO_UPA_corre(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,SNR,Kdir,K1,K2,f,D,dist_ris,lr,lt,Pt, no_iter,no_mat) 
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

sin_inc = dist_ris/sqrt( dist_ris^2+lt^2);
sin_ref = -(D-dist_ris)/sqrt( (D-dist_ris)^2+lr^2);

abnormal_phase = -k*(sin_inc+sin_ref)*RISPosition(1,:);


[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_UPA_rec_ob_corre(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,no_mat,Kdir,K1,K2,f,dist_ris,alpha_dir);

[Hdirt_K1000,H1t_K1000,H2t_K1000] = chan_mat_RIS_new_model_UPA_rec_ob_corre(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,1,100000,100000,100000, f,dist_ris,alpha_dir);
Hdir_K1000 = Hdirt_K1000{1}; 
H1_K1000 = H1t_K1000{1}; 
H2_K1000 = H2t_K1000{1};
Hdir1 = zeros(size(Hdir_K1000));
Qinit = eye(Nt)*(Pt/Nt);
omega_init = ones(1,Nris); 
c = 10;
[dC6,~,~,RIS_phase_opt_K1000] = GPM_Stefan(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1_K1000*sqrt(SNR),H2_K1000,500,Qinit*c^2,omega_init/c,c);
RIS_phase_K1000 =   angle(RIS_phase_opt_K1000);
 
offset = 0;
[Hdirt_siso,H1t_siso,H2t_siso] = chan_mat_RIS_new_model_UPA_rec_ob_corre(1,1,1,1,Nris_x,Nris_y,lt+lambda*offset,lr-lambda*offset,D,1,1000000,1000000,1000000,f,dist_ris,3);
Hdir_siso = Hdirt_siso{1}; H1_siso = H1t_siso{1}; H2_siso = H2t_siso{1};
Phase_RIS_focusing = -angle( H1_siso .* transpose( H2_siso ) );

AchievableRate_anom_hold = 0;
AchievableRate_K1000_solution_hold = 0;
AchievableRate_focusing_hold = 0;
AchievableRate_hold = 0;

phase_output = [ ];
SingleValue_output = [ ];

for i = 1:no_mat 
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
    
    StreamNum = min(Nt,Nr);
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); 
    c = 10;
    Hdir1 = zeros(size(Hdir));
    [dC6,~,~,RIS_phase_opt] = GPM_Stefan(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    AchievableRate_hold = AchievableRate_hold+dC6(end);
%     AchievableRate_hold = AchievableRate_hold+dC6;
    RIS_config_OptPhase = exp(1j*angle(RIS_phase_opt) );
    RIS_config_LoSOpt = exp(1j*RIS_phase_K1000);
    RIS_config_abnormal = exp(1j*abnormal_phase);
    RIS_config_focusing = exp(1j*Phase_RIS_focusing); 
end
