function rate_vs_time(Nt,Nr,Nris,SNR,K,f,D,dist_ris,lr,lt,Pt,nonDIR) 
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
no_mat = 200;
% Pt = 1; 
quant_bits = 0;
no_iter = 1000;
alpha_dir = 3;

[Hdirt,H1t,H2t] = chan_mat_RIS_surf_univ_new(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);

no_iter_pgm = no_iter+500;
C1 = 0;
C2 = zeros(1,no_iter+1);
C4 = zeros(1,no_iter_pgm+1);
L = 100;
time_pgm = zeros(1,no_iter_pgm+1);
time_ao = zeros(1,no_iter+1);

for i = 1:no_mat 
    Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
    
    if nonDIR
        Hdir = zeros(Nr,Nt);
        c = 10;
    else
        c = sqrt(norm(Hdir)/norm(H2*H1))*max(sqrt(Pt),1)/sqrt(Pt)*10; 
    end
    
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); % exp(-1i*rand(Nris,1)*2*pi);
    
    [dC1,dC2,omega_init1,Qinit1,iter_time_ao] = ...
        cap_zhang_mod(Nt,Nr,Nris,SNR,L,quant_bits,Pt,Hdir,H1,H2,no_iter);
    C1 = C1+dC1;
    C2 = C2+dC2;    
    [dC4,pgm_step(i,:),iter_time_pgm] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter_pgm,Qinit*c^2,omega_init/c,c);
    C4 = C4+dC4;
    time_pgm = time_pgm+iter_time_pgm;
    time_ao = time_ao+iter_time_ao;
end
Cao = C2/no_mat;
Cpgm = C4/no_mat;
time_ao = time_ao/no_mat*1000;
time_pgm = time_pgm/no_mat*1000;
filename = 'rate_vs_time_with_dir';
delete(filename);
save(filename,'Cao','time_ao','Cpgm','time_pgm');

plot(time_ao,Cao,'r-','DisplayName','AO'); hold on;
plot(time_pgm,Cpgm,'b-','DisplayName','PGM'); hold off;
xlabel('Time [ms]'); ylabel('Achievable rate [bit/s/Hz]');
xlim([0 500]);
legend('show','Location','SouthEast');

