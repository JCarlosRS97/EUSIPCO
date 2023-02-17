function [Cout,f_arr] = RateVsFreq(nonDIR)
Nt = 8;
Nr = 4;
SNR = db2pow(120);
K = 1;
% nonDIR = 0;
% D = 1200;
dist_ris = 40;
ris_size = 1;
D = 500; 
% nonDIR = 1;

lt = 20;
lr = 100;
% no_mat = 3;
Pt = 1; 
quant_bits = 0;
no_iter = 200;
alpha_dir = 3;

f_arr = [0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 9 10]*1e9;
% f_arr = [0.5 0.75 1]*1e9;
% f_arr = [10]*1e9;
% f = [0.5 0.75 1 1.5 2 2.5 3]*1e9;

if nonDIR==0
    filename = 'RateVsFreq_DIR';
else
    filename = 'RateVsFreq_nonDIR';
end
% filename = strcat(filename,'.mat'); 
    
lambda = 3e8./f_arr;
Nris_arr = floor(ris_size./(lambda/2)).^2;

for k = 1:length(f_arr)
    f = f_arr(k)
    Nris = Nris_arr(k);
    C = zeros(1,no_iter+1);
    Qinit = eye(Nt)*(Pt/Nt);
    omega_init = ones(1,Nris); 
    
    if f<=5e9
        no_mat = 10;
    elseif f<=7e9
        no_mat = 5;
    else
        no_mat = 2;
    end
    [Hdirt,H1t,H2t] = chan_mat_RIS_surf_univ_new(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,alpha_dir);
    for i = 1:no_mat
        Hdir = Hdirt{i}; H1 = H1t{i}; H2 = H2t{i};
        if nonDIR
            Hdir = zeros(Nr,Nt);
            c = 10;
        else
%             c = sqrt(norm(Hdir,'fro')/norm(H2*H1,'fro'))*sqrt(10);
            c = sqrt(norm(Hdir)/norm(H2*H1))*max(sqrt(Pt),1)/sqrt(Pt)*10;
        end
        [dC,~] = GPM_rescale2(Nt,Nr,Nris,1,quant_bits,Pt,Hdir*sqrt(SNR)/c,H1*sqrt(SNR),H2,no_iter,Qinit*c^2,omega_init/c,c);
        C = C+dC;
    end
    C = C/no_mat;
%     plot(C); hold on; 
    Cout(k) = max(C);
end
delete(filename);
save(filename,'f_arr','Cout');

% hold off;
% plot(f_arr/1e9,Cout);
% xlabel('$f$ [GHz]'); ylabel('Achievable rate [bpcu]');
% end
