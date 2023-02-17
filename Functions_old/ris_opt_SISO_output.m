function [ARoutput,AR_K1000_output,AR_anom_output, phase_output] = ris_opt_SISO_output(Nris,SNR,K,f,D,dist_ris,lr,lt,Pt) 
% clear; close all;
% Nt = 8;
% Nr = 4;
% Nris = 15^2;
% SNR = db2pow(120);
% K = 1;
Nt=1;
Nr=1;
no_mat = 30;
alpha_dir = 3;

lambda = 3e8/f;     % wavelength
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber
RISPosition = zeros(3,Nris);
N1 = sqrt(Nris);
for i = 1:N1
    for ii = 1:N1
        RISPosition( 1, (i-1)*N1+ii ) = (i - (N1+1)/2)*dris;
        RISPosition( 2, (i-1)*N1+ii ) = (ii - (N1+1)/2)*dris;
    end
end
ris_arr = RISPosition;
for l1 = 1:Nris
    d1_ex(l1) = dist_ris/sqrt(lt^2+dist_ris^2) * ris_arr(1,l1) ;
    d2_ex(l1) = -(D-dist_ris)/sqrt(lt^2+(D-dist_ris)^2) * ris_arr(1,l1) ;
end

abnormal_phase = -k*( d1_ex + d2_ex );


[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_ob(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

[Hdirt_K1000,H1t_K1000,H2t_K1000] = chan_mat_RIS_new_model_ob(Nt,Nr,Nris,lt,lr,D,1,1000,f,dist_ris,alpha_dir);
Hdir_K1000 = Hdirt_K1000{1}; 
H1_K1000 = H1t_K1000{1}; 
H2_K1000 = H2t_K1000{1};
RIS_phase_K1000 = -angle(diag(H2_K1000)*H1_K1000);
 
AchievableRate_anom_hold = 0;
AchievableRate_K1000_solution_hold = 0;
AchievableRate_hold = 0;

phase_output = [ ];
for i = 1:no_mat 
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
    
    RIS_phase = -angle(diag(H2)*H1);
    myomega = exp(1j*RIS_phase);
    Z = sqrt(SNR)*(H2*diag(myomega)*H1);
    AchievableRate = real(log(det(1+Z*Pt*Z')))/log(2);
    RIS_phase_opt = myomega;
    AchievableRate_hold = AchievableRate_hold+AchievableRate;
    phase_output = [ phase_output; transpose( RIS_phase_opt)];
    
    myomega = exp(1j*RIS_phase_K1000);
    Z = sqrt(SNR)*(H2*diag(myomega)*H1);
    AchievableRate_K1000_solution = real(log(det(1+Z*Pt*Z')))/log(2);
    
    AchievableRate_K1000_solution_hold = AchievableRate_K1000_solution_hold+AchievableRate_K1000_solution;
    
    RIS_phase_anom = abnormal_phase;
    myomega = exp(1j*RIS_phase_anom);
    Z = sqrt(SNR)*(H2*diag(myomega)*H1);
    AchievableRate_anom = real(log(det(1+Z*Pt*Z')))/log(2);
    AchievableRate_anom_hold = AchievableRate_anom_hold + AchievableRate_anom;
  
end

ARoutput = AchievableRate_hold/no_mat;
AR_K1000_output = AchievableRate_K1000_solution_hold/no_mat;
AR_anom_output = AchievableRate_anom_hold/no_mat;
