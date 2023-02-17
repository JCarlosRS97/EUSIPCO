function [ARoutput,AR_K1000_output,AR_anom_output, AR_focusing_output, phase_output, ...
    SingleValue_output, LeftSingleVector, RightSingleVector ] = ris_opt_MIMO_output_multi(Nt,Nr,Nris_x,Nris_y,SNR,K,f,D,dist_ris,lr,lt,Pt, no_iter,no_mat) 
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


[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_rec_ob(Nt,Nr,Nris_x,Nris_y,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

[Hdirt_K1000,H1t_K1000,H2t_K1000] = chan_mat_RIS_new_model_rec_ob(Nt,Nr,Nris_x,Nris_y,lt,lr,D,1,100000,f,dist_ris,alpha_dir);
Hdir_K1000 = Hdirt_K1000{1}; 
H1_K1000 = H1t_K1000{1}; 
H2_K1000 = H2t_K1000{1};
Hdir1 = zeros(size(Hdir_K1000));
Qinit = eye(Nt)*(Pt/Nt);
omega_init = ones(1,Nris); 
c = 10;
[dC6,~,~,RIS_phase_opt_K1000] = GPM_Stefan(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1_K1000*sqrt(SNR),H2_K1000,500,Qinit*c^2,omega_init/c,c);
RIS_phase_K1000 =   angle(RIS_phase_opt_K1000);
% dC6(end)
% RIS_phase_K1000 = -angle(diag(H2_K1000)*H1_K1000);
 
[Hdirt_siso,H1t_siso,H2t_siso] = chan_mat_RIS_new_model_rec_ob(1,1,Nris_x,Nris_y,lt,lr,D,1,1000000,f,dist_ris,3);
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
    

    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); 
    c = 10;
    Hdir1 = zeros(size(Hdir));
    [dC6,~,~,RIS_phase_opt] = GPM_Stefan(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    AchievableRate_hold = AchievableRate_hold+dC6(end);
%     AchievableRate_hold = AchievableRate_hold+dC6;
    phase_output = [ phase_output; transpose( RIS_phase_opt)];
    H_comp = Hdir1 + H2 * diag(exp(1j*angle(RIS_phase_opt) ) )* H1;
    [U,S,V] = svd(H_comp);
    SingleValue_output = [ SingleValue_output; (diag(S))'  ];
    LeftSingleVector{i} = U;
    RightSingleVector{i} = V;
    
    
    omega_init = exp(1j*RIS_phase_K1000);
    [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    AchievableRate_K1000_solution_hold = AchievableRate_K1000_solution_hold+dC6(end);

    RIS_phase_anom = abnormal_phase;
    omega_init = exp(1j*RIS_phase_anom);
    [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    
    AchievableRate_anom_hold = AchievableRate_anom_hold + dC6(end) ;
    
    
    omega_init = exp(1j*Phase_RIS_focusing);
    [dC6,~,~,RIS_phase_out] = GPM_FixedPhase(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,5,Qinit*c^2,omega_init/c,c);
    AchievableRate_focusing_hold = AchievableRate_focusing_hold+dC6(end);
  
end
phase_output = [ phase_output; transpose( RIS_phase_opt_K1000)];
phase_output = [ phase_output; ( exp(1j*RIS_phase_anom))];

ARoutput = AchievableRate_hold/no_mat;
AR_K1000_output = AchievableRate_K1000_solution_hold/no_mat;
AR_anom_output = AchievableRate_anom_hold/no_mat;
AR_focusing_output = AchievableRate_focusing_hold/ no_mat;