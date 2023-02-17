function plot_ris_cap_onlyindir(Nt,Nr,Nris,SNR,K,f,D,dist_ris,lr,lt,Pt) 
% clear; close all;
% Nt = 8;
% Nr = 4;
% Nris = 15^2;
% SNR = db2pow(120);
% K = 1;
nonDIR = 1;
% D = 500;
% dist_ris = D-40;
% f = 2000e6;

% lt = 20;
% lr = 100;
% h = 100;
% D = 50;
no_mat = 2;
% Pt = 1; 
quant_bits = 0;
no_iter = 500;
alpha_dir = 3;

[Hdirt,H1t,H2t] = chan_mat_RIS_new_model(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);


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
        c = 1/sqrt(norm(Hdir,'fro')/norm(H2*H1,'fro'))*sqrt(10);
    end
    
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); % exp(-1i*rand(Nris,1)*2*pi);
    
    [dC1,dC2,omega_init1,Qinit1] = ...
        cap_zhang_mod(Nt,Nr,Nris,SNR,L,quant_bits,Pt,Hdir,H1,H2,no_iter);
    C1 = C1+dC1;
    C2 = C2+dC2;    
    [dC4,pgm_step(i,:)] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    C4 = C4+dC4;
    
%     [dC5,apgm_step1(i,:),apgm_step2(i,:)] = GPM_acc_rescale_modified2(Nt,Nr,Nris,lt ,lr ,1,D,quant_bits,Pt,...
%          Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
%     C5 = C5+dC5; 
  
end
figure()
semilogx(1:length(C2),C2/no_mat,'g','DisplayName','AO'); hold on;
semilogx(1:length(C4),C4/no_mat,'r','DisplayName','PGM'); hold on;
% semilogx(1:length(C5),C5/no_mat,'b','DisplayName','APGM'); hold on;

xlabel('Iteration number'); ylabel('Achievable rate [bit/s/Hz]');
xlim([0 no_iter]); 
legend('show','Location','SouthEast');


