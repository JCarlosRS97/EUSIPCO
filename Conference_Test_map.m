% function dif_dist_BS_fig(fig)
addpath('Functions/')
addpath(genpath('Prol_1D-master/'))
clear
clc

Nt_y = 8;
Nt_z = 8;
Nr_y = 8;
Nr_z = 8;
Nris_x = 32;
Nris_y = 32;

N_0_dBm = -170;
BW = 20e6;

N_0 = 10.^((N_0_dBm - 30)/10);
P_noise = N_0*BW;
Pt = 10^((-20 - 30)/10);

f = 3.5e9;
lambda = 3e8/f;

r0_RIS = [0,0,0];
r0_TX = [5,-5,0];

x = [0.5:0.2:16];
y = [-10.1:0.2:10.1];

[X,Y] = meshgrid(x,y);
X_vec = X(:);
Y_vec = Y(:);

%%

output_ND_num_WA = zeros(size(X_vec));
output_FOC_num_WA = zeros(size(X_vec));
output_ND_PSWF_WA = zeros(size(X_vec));
output_FOC_PSWF_WA = zeros(size(X_vec));
output_ND_num_EQ = zeros(size(X_vec));
output_FOC_num_EQ = zeros(size(X_vec));
output_ND_PSWF_EQ = zeros(size(X_vec));
output_FOC_PSWF_EQ = zeros(size(X_vec));

for indexN = 1:length(X_vec)
    %     d_RX_n = (D_dist( indexN ) );
    %
    %     d_ris = d_RX_n/2;


    r0_RX = [X_vec(indexN),Y_vec(indexN),0];

    [rate_ND_NUM_WA, rate_ND_PSWF_WA, rate_FOC_NUM_WA, rate_FOC_PSWF_WA,...
        rate_ND_NUM_EQ, rate_ND_PSWF_EQ, rate_FOC_NUM_EQ, rate_FOC_PSWF_EQ] = ...
        ris_opt_MIMO_UPA_aux(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,Pt,P_noise,f,r0_TX, r0_RX,r0_RIS);

    output_ND_num_WA(indexN) = rate_ND_NUM_WA;
    output_ND_PSWF_WA(indexN) = rate_ND_PSWF_WA;
    output_FOC_num_WA(indexN) = rate_FOC_NUM_WA;
    output_FOC_PSWF_WA(indexN) = rate_FOC_PSWF_WA;
    output_ND_num_EQ(indexN) = rate_ND_NUM_EQ;
    output_ND_PSWF_EQ(indexN) = rate_ND_PSWF_EQ;
    output_FOC_num_EQ(indexN) = rate_FOC_NUM_EQ;
    output_FOC_PSWF_EQ(indexN) = rate_FOC_PSWF_EQ;



end
%%
filename = [ 'Results\Map_RIS_MIMO_Nty' num2str(Nt_y) 'Ntz' num2str(Nt_z) 'Nry' num2str(Nr_y) 'Nrz' num2str(Nr_z)...
    'Nrisx' num2str(Nris_x)  '_drx.mat'];
save( filename )
%%
close all
rate_ND_PSWF = output_ND_PSWF_WA./output_ND_num_WA;
rate_FOC_num = output_FOC_num_WA./output_ND_num_WA;
rate_FOC_PSWF = output_FOC_PSWF_WA./output_ND_num_WA;

figure(1),surf(X,Y, reshape(rate_ND_PSWF, [21,16])),view(2),caxis([0.9 1]);
figure(2),surf(X,Y, reshape(rate_FOC_num, [21,16])),view(2),caxis([0.9 1]);
figure(3),surf(X,Y, reshape(rate_FOC_PSWF, [21,16])),view(2),caxis([0.9 1]);

[f_ND_PSWF,x_ND_PSWF] = ecdf(output_ND_PSWF);
[f_ND_num,x_ND_num] = ecdf(output_ND_num);
[f_FOC_num,x_FOC_num] = ecdf(output_FOC_num);
[f_FOC_PSWF,x_FOC_PSWF] = ecdf(output_FOC_PSWF);

figure()
plot(x_ND_num,1-f_ND_num),hold on
plot(x_ND_PSWF,1-f_ND_PSWF)
plot(x_FOC_num,1-f_FOC_num)
plot(x_FOC_PSWF,1-f_FOC_PSWF)
legend('ND-num','ND-PWSF','FOC-num','FOC-PSWF')