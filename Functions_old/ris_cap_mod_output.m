function [output, phase_output] = ris_cap_mod_output(Nt,Nr,Nris,SNR,K,f,D,dist_ris,lr,lt,Pt) 
% clear; close all;
% Nt = 8;
% Nr = 4;
% Nris = 15^2;
% SNR = db2pow(120);
% K = 1;
nonDIR = 0;
no_mat = 30;
quant_bits = 0;
no_iter = 100;
alpha_dir = 3;

[Hdirt,H1t,H2t] = chan_mat_RIS_new_model_ob(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);


C1 = 0;
C2 = zeros(1,no_iter+1);
C3 = zeros(1,no_iter+1);
C4 = zeros(1,no_iter+1);
C5 = zeros(1,no_iter+1);
C6 = zeros(1,no_iter+1);
C7 = zeros(1,no_iter+1);
L = 100;

for i = 1:no_mat 
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
    
    if nonDIR
        Hdir = zeros(Nr,Nt);
        c = 10;
    else
        %         c = sqrt(norm(Hdir)/norm(H2*H1));
%         c = 1/sqrt(norm(Hdir,'fro')/norm(H2*H1,'fro'))*sqrt(10);
        c = sqrt(norm(Hdir)/norm(H2*H1))*max(sqrt(Pt),1)/sqrt(Pt)*10; 
    end
    
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); % exp(-1i*rand(Nris,1)*2*pi);
    
%     [dC1,dC2,omega_init1,Qinit1] = ...
%         cap_zhang_mod(Nt,Nr,Nris,SNR,L,quant_bits,Pt,Hdir,H1,H2,no_iter);
%     C1 = C1+dC1;
%     C2 = C2+dC2;    
    
%     [dC4,pgm_step(i,:)] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
%     C4 = C4+dC4;
    
%     [dC5,apgm_step1(i,:),apgm_step2(i,:)] = GPM_acc_rescale_modified2(Nt,Nr,Nris,lt ,lr ,1,D,quant_bits,Pt,...
%          Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
%     C5 = C5+dC5; 
    % only INDIR
    c = 10;
    Hdir1 = zeros(size(Hdir));
%     [dC6,~,~] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    [dC6,~,~,phase_output] = GPM_Stefan(Nt,Nr,Nris,1,quant_bits,Pt,Hdir1*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);

    C6 = C6+dC6;  
    % only  DIR
%     c=1;
%     H11= zeros(size(H1));
%     H21= zeros(size(H2));     
%     [dC7] = GPM_acc_Q(Nt,Nr,Nris,lt,lr,1,D,quant_bits,Pt,...
%          Hdir*sqrt(SNR)/c,H11*sqrt(SNR),H21,no_iter,Qinit*c^2,omega_init/c,c);
%     C7 = C7+dC7;     
    
end

output = C6/no_mat;
% figure()
% % C4 = C4/no_mat; C5 = C5/no_mat;
% % plot(1:length(C1),C1/no_mat,'b-s','DisplayName','C1'); hold on;
% semilogx(1:length(C2),C2/no_mat,'g','DisplayName','AO'); hold on;
% semilogx(1:length(C4),C4/no_mat,'r','DisplayName','PGM'); hold on;
% % semilogx(1:length(C5),C5/no_mat,'b','DisplayName','APGM'); hold on;
% semilogx(1:length(C6),C6/no_mat,'b--','DisplayName','Ind. link'); hold on;
% semilogx(1:length(C7),C7/no_mat,'r:','DisplayName','Dir. link'); hold on;
% xlabel('Iteration number'); ylabel('Achievable rate [bit/s/Hz]');
% xlim([0 no_iter]); 
% legend('show','Location','SouthEast');

% number of gradient computations (Ic), and additional number of updates and
% projections (Ip)
% PGM
% Ic_pgm = grad_no(C4(2:end))
% Ip_pgm = Ic_pgm+update_proj_no(C4(2:end),pgm_step)
% % APGM
% Ic_apgm = grad_no(C5(2:end))
% Ip1_apgm = Ic_apgm+update_proj_no(C5(2:end),apgm_step1)
% Ip2_apgm = Ic_apgm+update_proj_no(C5(2:end),apgm_step2)
% % end
