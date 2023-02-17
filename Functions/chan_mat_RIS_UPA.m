function [Hdir,H1,H2] = chan_mat_RIS_UPA(tx_arr, rx_arr, ris_arr,Kdir,K1,K2,f, no_mat)

lambda = 3e8/f;     % wavelength
k = 2*pi/lambda;    % wavenumber


Gt = 1;
Gr = 1;
G = 1;
Nt = size(tx_arr,2);
Nr = size(rx_arr,2);
Nris = size(ris_arr,2);

%% H1 LOS
% Constant1 = sqrt( dris*dris*lambda*lambda/64/pi^3) ;
Constant1 = sqrt( lambda^4/(256*pi^2)) ; 
H1_los = zeros(Nris,Nt);
for i = 1:Nt
    for ii = 1:Nris
        distance_nt_nm = norm( tx_arr(:,i) - ris_arr(:,ii) );
        cos_theta_t_nt_nm = tx_arr(3,i)/distance_nt_nm;
        Amplitude_nt_nm = Constant1*G*Gt*cos_theta_t_nt_nm/distance_nt_nm/distance_nt_nm;
        Phase_term_nt_nm = exp( 1j*k*distance_nt_nm);
%         H1_los(ii,i) = sqrt(Amplitude_nt_nm)*Phase_term_nt_nm;
        H1_los(ii,i) = Phase_term_nt_nm/distance_nt_nm/4/pi;

    end
end

%% H2 LOS
H2_los = zeros(Nr,Nris);
for i = 1:Nr
    for ii = 1:Nris
        distance_nr_nm = norm( rx_arr(:,i) - ris_arr(:,ii) );
        cos_theta_r_nr_nm = rx_arr(3,i)/distance_nr_nm;       
        Amplitude_nr_nm = Constant1*G*Gr*cos_theta_r_nr_nm/distance_nr_nm/distance_nr_nm;
        Phase_term_nr_nm = exp( 1j*k*distance_nr_nm);
%         H2_los(i,ii) = sqrt(Amplitude_nr_nm)*Phase_term_nr_nm;
        H2_los(i,ii) = Phase_term_nr_nm/distance_nr_nm/4/pi;
    end
end


for i1 = 1:Nr
    for j1 = 1:Nt
        d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
    end 
end

% Hdir_los = exp(1i*k*d) .* sqrt( Gt*Gr*(lambda/(4*pi))^2./(d.^(2)) );
Hdir_los = exp(1i*k*d) .* sqrt( Gt*Gr*(1/(4*pi))^2./(d.^(2)) );


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
