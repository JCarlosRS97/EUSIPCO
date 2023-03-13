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

d_0 = 8;
gamma = [5:5:85]*pi/180;

r0_RIS = [0,0,0];

%%

output_ND_num_WA = zeros(size(gamma));
output_FOC_num_WA = zeros(size(gamma));
output_ND_PSWF_WA = zeros(size(gamma));
output_FOC_PSWF_WA = zeros(size(gamma));
output_ND_num_EQ = zeros(size(gamma));
output_FOC_num_EQ = zeros(size(gamma));
output_ND_PSWF_EQ = zeros(size(gamma));
output_FOC_PSWF_EQ = zeros(size(gamma));

for indexN = 1:length(gamma)
    r0_TX = [d_0*sin(gamma(indexN)),-d_0*cos(gamma(indexN)),0];

    r0_RX = [d_0*sin(gamma(indexN)),d_0*cos(gamma(indexN)),0];

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
filename = [ 'Results\Angles_RIS_MIMO_Nty' num2str(Nt_y) 'Ntz' num2str(Nt_z) 'Nry' num2str(Nr_y) 'Nrz' num2str(Nr_z)...
    'Nrisx' num2str(Nris_x) 'd_0' num2str(d_0)  '.mat'];
save( filename )
%%
% close all

figure()
plot(gamma*180/pi,output_ND_num_WA),hold on
plot(gamma*180/pi,output_ND_PSWF_WA)
plot(gamma*180/pi,output_FOC_num_WA)
plot(gamma*180/pi,output_FOC_PSWF_WA)

plot(gamma*180/pi,output_ND_num_EQ),hold on
plot(gamma*180/pi,output_ND_PSWF_EQ)
plot(gamma*180/pi,output_FOC_num_EQ)
plot(gamma*180/pi,output_FOC_PSWF_EQ)
legend('ND-num','ND-PWSF','FOC-num','FOC-PSWF')