% function dif_dist_BS_fig(fig)
addpath('Functions/')
addpath(genpath('Prol_1D-master/'))
clear
clc

Nt_y = 20;
Nt_z = 20;
Nr_y = 20;
Nr_z_vec = 10:5:40;

N_0_dBm = -160;
BW = 20e6;

N_0 = 10.^((N_0_dBm - 30)/10);
P_noise = N_0*BW;
Pt = 10^((-20 - 30)/10);


f = 3.5e9;
lambda = 3e8/f;

r0_RIS = [0,0,0];
d_TX = 20;
d_RX = 20;

r0_TX = [-d_TX,0,0];
    r0_RX = [d_RX,0,0];

n_iter = 200;

%%
Indir_output_K100000_hold = [  ];
Indir_output_K100000_hold_K1000Solution = [  ];
Indir_output_K100000_hold_anom = [  ];
Indir_output_K100000_hold_focusing = [  ];

Nris_z = 20;
Nris_y = 20;


output_optimal = zeros(size(d_RX));
output_focusing = zeros(size(d_RX));
output_miller = zeros(size(d_RX));
output_approx = zeros(size(d_RX));

for indexN = 1:length(Nr_z_vec)
    Nr_z = (Nr_z_vec( indexN ) );

    [optimal_rate, focusing_rate, opt_rate, miller_rate, approx_rate] = ...
        ris_opt_MIMO_UPA_aux(Nt_y,Nt_z,Nr_y,Nr_z,Nris_z,Nris_y,P_noise,f,r0_TX, r0_RX, r0_RIS,Pt, n_iter);

    output_optimal(indexN) = optimal_rate;
    output_focusing(indexN) = focusing_rate;
    output_miller(indexN) = miller_rate;
    output_approx(indexN) = approx_rate;


end
%%
filename = [ 'Results\Rate_RIS_MIMO_Nty' num2str(Nt_y) 'Ntz' num2str(Nt_z) 'Nry' num2str(Nr_y) 'Nrz' num2str(Nr_z)...
    'Nrisy' num2str(Nris_z)  '_RIS.mat'];
save( filename )
%%
% surf(real(reshape(ObAtRIS_OptPhase_K100000{1}, [Nris_y,NRIS_oneside])))
figure(),plot(Nr_z*Nr_y, output_optimal), hold on
plot(Nr_z*Nr_y, output_focusing)
plot(Nr_z*Nr_y, output_miller)
plot(Nr_z*Nr_y, output_approx)
xlabel('Number of RIS elements'), ylabel('Spectral efficiency (bps/Hz)')
legend('Optimal scheme', 'Focusing scheme', 'Transparent scheme', 'Proposed scheme')