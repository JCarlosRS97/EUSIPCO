function [rate_ND_NUM_WA, rate_ND_PSWF_WA, rate_FOC_NUM_WA, rate_FOC_PSWF_WA,...
    rate_ND_NUM_EQ, rate_ND_PSWF_EQ, rate_FOC_NUM_EQ, rate_FOC_PSWF_EQ] =...
    ris_opt_MIMO_UPA_aux(Nt_x,Nt_z,Nr_x,Nr_z,Nris_y,Nris_z,Pt,P_noise,f,r0_TX, r0_RX,r0_RIS)

lambda = 3e8/f;     % wavelength
d_e = lambda/2;
k_0 = 2*pi/lambda;    % wavenumber
Nt = Nt_x * Nt_z;
Nr = Nr_x * Nr_z;
Nris = Nris_y*Nris_z;



%% Place antennas
% RIS array
ris_arr = zeros(3,Nris);
% [ris_v_x, ris_v_y]  = place_antennas(Nris_x, Nris_y, d_e);
% ris_arr(1,:) = ris_v_x + r0_RIS(1);
% ris_arr(2,:) = ris_v_y + r0_RIS(2);
% ris_arr(3,:) = r0_RIS(3);

[ris_v_y, ris_v_z]  = place_antennas(Nris_y, Nris_z, d_e);
ris_arr(1,:) = r0_RIS(1);
ris_arr(2,:) = ris_v_y + r0_RIS(2);
ris_arr(3,:) = ris_v_z + r0_RIS(3);

%TX array
tx_arr = zeros(3,Nt);
[tx_v_x, tx_v_z]  = place_antennas(Nt_x, Nt_z, d_e);
tx_arr(1,:) = tx_v_x + r0_TX(1);
tx_arr(2,:) = r0_TX(2);
tx_arr(3,:) = tx_v_z + r0_TX(3);

% RX array
rx_arr = zeros(3,Nr);
[rx_v_x, rx_v_z]  = place_antennas(Nr_x, Nr_z, d_e);
rx_arr(1,:) = rx_v_x + r0_RX(1);
rx_arr(2,:) = r0_RX(2);
rx_arr(3,:) = rx_v_z + r0_RX(3);

[~,H1,H2] = chan_mat_RIS_UPA(tx_arr, rx_arr, ris_arr,f);
Hdir1 = zeros(Nt,Nr);


%% Analytical computation

Delta_x_T = d_e*Nt_x/2;
Delta_z_T = d_e*Nt_z/2;
Delta_x_R = d_e*Nr_x/2;
Delta_z_R = d_e*Nr_z/2;
Delta_y_ris = d_e*Nris_y/2;
Delta_z_ris = d_e*Nris_z/2;
d_1 = norm(r0_RIS-r0_TX);
d_2 = norm(r0_RIS-r0_RX);

factor_t_e = abs(r0_TX(1))*abs(r0_TX(2))/d_1^2;
factor_r_e = abs(r0_RX(1))*abs(r0_RX(2))/d_2^2;

N_DoF_1 =  ceil((2*Delta_x_T)*(2*Delta_z_T)*(2*Delta_y_ris)*(2*Delta_z_ris)/(lambda^2*d_1^2)*factor_t_e);
N_DoF_2 =  ceil((2*Delta_x_R)*(2*Delta_z_R)*(2*Delta_y_ris)*(2*Delta_z_ris)/(lambda^2*d_2^2)*factor_r_e);
N_DoF = min(N_DoF_1,N_DoF_2);

c_t_e_x = k_0*Delta_x_T*Delta_y_ris/d_1*factor_t_e;
c_t_e_z = k_0*Delta_z_T*Delta_z_ris/d_1;

c_r_e_x = k_0*Delta_x_R*Delta_y_ris/d_2*factor_r_e;
c_r_e_z = k_0*Delta_z_R*Delta_z_ris/d_2;

tx_x_vec = unique(tx_v_x).'/Delta_x_T;
tx_z_vec = unique(tx_v_z).'/Delta_z_T;

rx_x_vec = unique(rx_v_x).'/Delta_x_R;
rx_z_vec = unique(rx_v_z).'/Delta_z_R;

ris_y_vec = unique(ris_v_y).'/Delta_y_ris;
ris_z_vec = unique(ris_v_z).'/Delta_z_ris;

% Patterns
[P_TX_RIS, eigenval_ord_1] = compute_analytical_pattern(c_t_e_x, c_t_e_z,...
    tx_x_vec,tx_z_vec, k_0,50);

[P_RIS_TX, ~] = compute_analytical_pattern(c_t_e_x, c_t_e_z,...
    ris_y_vec,ris_z_vec, k_0,50);
[P_RIS_RX, eigenval_ord_2] = compute_analytical_pattern(c_r_e_x, c_r_e_z,...
    ris_y_vec,ris_z_vec, k_0,50);
% pswfs_ris_in = compute_analytical_pattern(c_tx_x, c_tx_y, x,y);
% pswfs_ris_out = compute_analytical_pattern(c_rx_x, c_rx_y, x,y);
% [pswfs_rx,~] = compute_analytical_pattern(c_r_e_x, c_r_e_z,...
%     rx_x_vec,rx_z_vec, k_0,N_DoF_2+50);

% Focusing function
F_RIS_TX = compute_focusing_function(ris_arr, r0_TX, k_0);
F_RIS_RX = conj(compute_focusing_function(ris_arr, r0_RX, k_0));
F_TX_RIS = conj(compute_focusing_function(tx_arr, r0_RIS, k_0));
% F_RX_RIS = compute_focusing_function(rx_arr, r0_RIS, k_0);

U_H = diag(F_RIS_TX)*P_RIS_TX;
V_G = diag(F_RIS_RX)*P_RIS_RX;

V_H = diag(F_TX_RIS)*P_TX_RIS;
g_1_cont = (2*k_0)^-2*eigenval_ord_1.'/factor_t_e;
g_2_cont = (2*k_0)^-2*eigenval_ord_2.'/factor_r_e;

g_1 = g_1_cont/(lambda/2)^4;
g_2 = g_2_cont/(lambda/2)^4;

S_approx = zeros(size(V_H,2),1);
max_N_eigenval = min([length(g_1),length(g_2), size(V_H,2)]);
S_approx(1:max_N_eigenval) = g_1(1:max_N_eigenval).*g_2(1:max_N_eigenval);
% S_approx(1:size(V_H,2)) = g_1(1:size(V_H,2)).*g_2(1:size(V_H,2));


%% Non diagonal

[U1,S1,V1] = svd(H1);
[U2,S2,V2] = svd(H2);

Phi_ND_NUM = V2*U1';
H_comp = H2*Phi_ND_NUM*H1;

eigenval = diag(S1).^2.*diag(S2).^2;

% Watefilling
D_ND_NUM_WA = water_fill(Pt, eigenval/P_noise);
Q_ND_NUM_WA = V1*diag(D_ND_NUM_WA)*V1';
rate_ND_NUM_WA = real(log2(det(eye(Nr)+H_comp*Q_ND_NUM_WA*H_comp'/P_noise)));

% Equipower
D_ND_NUM_EQ = zeros(size(eigenval));
D_ND_NUM_EQ(1:N_DoF) = Pt/N_DoF;
Q_ND_NUM_EQ = V1*diag(D_ND_NUM_EQ)*V1';
rate_ND_NUM_EQ = real(log2(det(eye(Nr)+H_comp*Q_ND_NUM_EQ*H_comp'/P_noise)));


%% Non diagonal analytical configuration
% Ris configuration
Phi_ND_PSWF = V_G*U_H';
H_comp = H2*Phi_ND_PSWF*H1;


% Watefilling
D_PSWF_WA = water_fill(Pt, S_approx/P_noise);
Q_PSWF_WA = V_H*diag(D_PSWF_WA)*V_H';

rate_ND_PSWF_WA = real(log2(det(eye(Nr)+H_comp*Q_PSWF_WA*H_comp'/P_noise)));

% Equipower
D_PSWF_EQ = zeros(length(S_approx),1);
D_PSWF_EQ(1:N_DoF) = Pt/N_DoF;
Q_PSWF_EQ = V_H*diag(D_PSWF_EQ)*V_H';
rate_ND_PSWF_EQ = real(log2(det(eye(Nr)+H_comp*Q_PSWF_EQ*H_comp'/P_noise)));

%% Focusing configurations

Phi_FOC_NUM = diag(F_RIS_RX)*diag(F_RIS_TX)';
H_comp = H2*Phi_FOC_NUM*H1;

[U,S,V] = svd(H_comp);

% Watefilling
D_FOC_NUM_WA = water_fill(Pt, diag(S).^2/P_noise);
Q_FOC_NUM_WA = V*diag(D_FOC_NUM_WA)*V';
rate_FOC_NUM_WA = real(log2(det(eye(Nr) + H_comp*Q_FOC_NUM_WA*H_comp'/P_noise)));
rate_FOC_PSWF_WA = real(log2(det(eye(Nr) + H_comp*Q_PSWF_WA*H_comp'/P_noise)));

% Equipower
D_FOC_NUM_EQ = zeros(length(S),1);
D_FOC_NUM_EQ(1:N_DoF) = Pt/N_DoF;
Q_FOC_NUM_EQ = V*diag(D_FOC_NUM_EQ)*V';
rate_FOC_NUM_EQ = real(log2(det(eye(Nr) + H_comp*Q_FOC_NUM_EQ*H_comp'/P_noise)));
rate_FOC_PSWF_EQ = real(log2(det(eye(Nr) + H_comp*Q_PSWF_EQ*H_comp'/P_noise)));

%% Opt configuration
% Qinit = eye(Nt)*(Pt/Nt);
% omega_init = ones(1,Nris);
% c=1;
% [dC6,stepsize1,stepsize2,~,RIS_phase_opt] = GPM_Stefan_carlos(Nt,Nr,Nris,1,0,Pt,Hdir1*sqrt(1/P_noise)/c,H1*sqrt(1/P_noise),H2,n_iter,Qinit*c^2,omega_init/c,c);
% 
% H_comp = H2 * diag(exp(1j*angle(RIS_phase_opt) ) )* H1;
% % [U,S,V] = svd(H_comp);
% 
% rate_D_OPT = dC6(end);

%% Miller configuration
% Delta_x_T = d_e*(Nt_x - 1)/2;
% Delta_z_T = d_e*(Nt_z - 1)/2;
% Delta_x_R = d_e*(Nr_x - 1)/2;
% Delta_z_R = d_e*(Nr_z - 1)/2;
% Delta_y_ris = d_e*(Nris_y - 1)/2;
% Delta_z_ris = d_e*(Nris_z - 1)/2;
% d_0 = norm(r0_RX-r0_TX);
% 
% c_y = k_0*Delta_x_T*Delta_x_R/d_0;
% c_z = k_0*Delta_z_T*Delta_z_R/d_0;
% 
% tx_y_vec = unique(tx_arr(2,:)).'/Delta_x_T;
% tx_z_vec = unique(tx_arr(3,:)).'/Delta_z_T;
% 
% rx_y_vec = unique(rx_arr(2,:)).'/Delta_x_R;
% rx_z_vec = unique(rx_arr(3,:)).'/Delta_z_R;
% % Patterns
% [pswfs_tx, g_2_cont] = compute_analytical_pattern(c_y, c_z,...
%     tx_y_vec,tx_z_vec, k_0);
% % pswfs_ris_in = compute_analytical_pattern(c_tx_x, c_tx_y, x,y);
% % pswfs_ris_out = compute_analytical_pattern(c_rx_x, c_rx_y, x,y);
% % pswfs_rx = compute_analytical_pattern(c_tx_x, c_tx_y,...
% %     rx_y_vec,rx_z_vec);
% 
% % Focusing function
% F_Tx = compute_focusing_function(tx_arr, r0_RX, k_0);
% F_Rx = compute_focusing_function(rx_arr, r0_TX, k_0);
% 
% input_pattern = F_Tx.*pswfs_tx;
% g_2 = g_2_cont/(lambda/2)^4;
% 
% S_miller = zeros(size(input_pattern,1),1);
% S_miller(1:length(g_2)) = g_2;
% % Watefilling
% D_miller = water_fill(Pt, S_miller/P_noise);
% Q_miller = input_pattern*diag(D_miller)*input_pattern';
% 
% % Ris configuration
% RIS_Phi_FOC = diag(conj(F_R1.*F_R2));
% H_comp = H2*RIS_Phi_FOC*H1;
% 
% miller_rate = real(log2(det(eye(Nr)+H_comp*Q_miller*H_comp'/P_noise)));


