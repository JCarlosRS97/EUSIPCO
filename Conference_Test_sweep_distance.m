% function dif_dist_BS_fig(fig)
addpath('Functions/')
addpath(genpath('Prol_1D-master/'))
clear
clc

Nt_y = 30;
Nt_z = 30;
Nr_y = 30;
Nr_z = 30;
Nris_x = 30;
Nris_y = 30;

N_0_dBm = -160;
BW = 20e6;

N_0 = 10.^((N_0_dBm - 30)/10);
P_noise = N_0*BW;
Pt = 10^((-20 - 30)/10);

f = 3.5e9;
lambda = 3e8/f;

d_ris = 10;
D_dist = 20;
lt = 10;
lr = 10;

n_iter = 200;

%%

output_optimal = zeros(size(D_dist));
output_focusing = zeros(size(D_dist));
output_miller = zeros(size(D_dist));
output_approx = zeros(size(D_dist));

for indexN = 1:length(D_dist)
    d_RX_n = (D_dist( indexN ) );

    [optimal_rate, focusing_rate, opt_rate, miller_rate, approx_rate] = ...
        ris_opt_MIMO_UPA_aux(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,P_noise,f,lt, lr, D_dist,d_ris,Pt, n_iter);

    output_optimal(indexN) = optimal_rate;
    output_focusing(indexN) = focusing_rate;
    output_miller(indexN) = miller_rate;
    output_approx(indexN) = approx_rate;


end
%%
filename = [ 'Results\Modes_RIS_MIMO_Nty' num2str(Nt_y) 'Ntz' num2str(Nt_z) 'Nry' num2str(Nr_y) 'Nrz' num2str(Nr_z)...
    'Nrisx' num2str(Nris_x)  '_drx.mat'];
save( filename )
%%
% surf(real(reshape(ObAtRIS_OptPhase_K100000{1}, [Nris_y,NRIS_oneside])))
figure(),plot(d_RX, output_optimal), hold on
plot(d_RX, output_focusing)
plot(d_RX, output_miller)
plot(d_RX, output_approx)
xlabel('$d_r$ (m)', 'Interpreter','latex','FontSize',10),
ylabel('Spectral efficiency (bps/Hz)', 'Interpreter','latex','FontSize',10)
legend('Optimal scheme', 'Focusing scheme', 'Transparent scheme',...
    'Proposed scheme', 'Interpreter','latex','FontSize',10)
axis([5,60,30,110])