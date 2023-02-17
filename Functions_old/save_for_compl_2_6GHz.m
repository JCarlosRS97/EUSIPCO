function save_for_compl_2_6GHz(Nt,Nr,Nris,SNR,K,f,D,dist_ris,lr,lt,Pt,nonDIR,no_iter) 
% clear; close all;
% Nt = 8;
% Nr = 4;
% Nris = 15^2;
% SNR = db2pow(120);
% K = 1;
% nonDIR = 0;
% D = 500;
% dist_ris = D-40;
% f = 2000e6;

% lt = 20;
% lr = 100;
% h = 100;
% D = 50;
no_mat = 20;
% Pt = 1; 
quant_bits = 0;
% no_iter = 500;
alpha_dir = 3;

[Hdirt,H1t,H2t] = chan_mat_RIS_surf_univ_new_pan(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);


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
        c = sqrt(norm(Hdir)/norm(H2*H1))*max(sqrt(Pt),1)/sqrt(Pt)*10; 
    end
    
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); 
    
    [~,dC2,omega_init1,Qinit1] = ...
        cap_zhang_mod(Nt,Nr,Nris,SNR,L,quant_bits,Pt,Hdir,H1,H2,no_iter);
    C2 = C2+dC2;    
    [dC4,pgm_step(i,:)] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
    C4 = C4+dC4;
       
    
end

Cao = C2/no_mat;
pgm_iter = mean(pgm_step,1);
Cpgm = C4/no_mat;

if nonDIR
    filename = sprintf('nonDIR_Nris_%d',Nris);
else
    filename = sprintf('DIR_Nris_%d',Nris);
end

delete(filename);
save(filename,'Cpgm','pgm_iter','Cao');
end
