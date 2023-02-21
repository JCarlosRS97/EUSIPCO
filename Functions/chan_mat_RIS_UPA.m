function [Hdir,H1,H2] = chan_mat_RIS_UPA(tx_arr, rx_arr, ris_arr,f)

lambda = 3e8/f;     % wavelength
k_0 = 2*pi/lambda;    % wavenumber


Gt = 1;
Gr = 1;
G = 1;
Nt = size(tx_arr,2);
Nr = size(rx_arr,2);
Nris = size(ris_arr,2);

%% H1 LOS
d_1 = sqrt((ris_arr(1,:).'-tx_arr(1,:)).^2 + (ris_arr(2,:).'-tx_arr(2,:)).^2 + (ris_arr(3,:).'-tx_arr(3,:)).^2);
H1_los = exp(1j*k_0*d_1)./(4*pi*d_1);

% Constant1 = sqrt( lambda^4/(256*pi^2)) ; 
% H1_los = zeros(Nris,Nt);
% for i = 1:Nt
%     for ii = 1:Nris
%         distance_nt_nm = norm( tx_arr(:,i) - ris_arr(:,ii) );
%         cos_theta_t_nt_nm = tx_arr(3,i)/distance_nt_nm;
% %         Amplitude_nt_nm = Constant1*G*Gt*cos_theta_t_nt_nm/distance_nt_nm/distance_nt_nm;
%         Phase_term_nt_nm = exp( 1j*k*distance_nt_nm);
% %         H1_los(ii,i) = sqrt(Amplitude_nt_nm)*Phase_term_nt_nm;
%         H1_los(ii,i) = Phase_term_nt_nm/distance_nt_nm/4/pi;
% 
%     end
% end



%% H2 LOS

d_2 = sqrt((rx_arr(1,:).'-ris_arr(1,:)).^2 + (rx_arr(2,:).'-ris_arr(2,:)).^2 + (rx_arr(3,:).'-ris_arr(3,:)).^2);
H2_los = exp(1j*k_0*d_2)./(4*pi*d_2);

% H2_los = zeros(Nr,Nris);
% for i = 1:Nr
%     for ii = 1:Nris
%         distance_nr_nm = norm( rx_arr(:,i) - ris_arr(:,ii) );
%         cos_theta_r_nr_nm = rx_arr(3,i)/distance_nr_nm;       
% %         Amplitude_nr_nm = Constant1*G*Gr*cos_theta_r_nr_nm/distance_nr_nm/distance_nr_nm;
%         Phase_term_nr_nm = exp( 1j*k*distance_nr_nm);
% %         H2_los(i,ii) = sqrt(Amplitude_nr_nm)*Phase_term_nr_nm;
%         H2_los(i,ii) = Phase_term_nr_nm/distance_nr_nm/4/pi;
%     end
% end

d = sqrt((rx_arr(1,:).'-tx_arr(1,:)).^2 + (rx_arr(2,:).'-tx_arr(2,:)).^2 + (rx_arr(3,:).'-tx_arr(3,:)).^2);
Hdir_los = exp(1j*k_0*d)./(4*pi*d);


% Hdir_los = exp(1i*k*d) .* sqrt( Gt*Gr*(lambda/(4*pi))^2./(d.^(2)) );
% Hdir_los = exp(1i*k*d) .* sqrt( Gt*Gr*(1/(4*pi))^2./(d.^(2)) );


Hdir = Hdir_los;
H1 = H1_los;
H2 = H2_los;
% Hdir = Rician_ewise(Hdir_los,no_mat,Kdir);
% H1 = Rician_ewise(H1_los,no_mat,K1);
% H2 = Rician_ewise(H2_los,no_mat,K2);
% 
% function Hout = Rician_ewise(Hlos,no_mat,K)
%     for ind = 1:no_mat
%         Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
%         Htot = (Hlos*sqrt(K)+Hnlos.*Hlos)/sqrt(K+1);
%         Hout{ind} = Htot;
%     end
% end

end
