function [Hdir,H1,H2] = chan_mat_RIS_new_model(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,varargin)

lambda = 3e8/f;     % wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber


Gt = 2;
Gr = 2;
G = 1;


% geometrical placement
% x, y and z axis 
% TX array
tx_arr = zeros(3,Nt); 
tx_arr(1,:) = -dist_ris;
tx_arr(3,:) = (1:Nt)*dt-(Nt+1)*dt/2+lt;
% RX array
rx_arr = zeros(3,Nr); 
rx_arr(1,:) = D-dist_ris;
rx_arr(3,:) = (1:Nr)*dr-(Nr+1)*dr/2+lr; 
% RIS
RISPosition = zeros(3,Nris);
N1 = sqrt(Nris);
for i = 1:N1
    for ii = 1:N1
        RISPosition( 1, (i-1)*N1+ii ) = (i - (N1+1)/2)*dris;
        RISPosition( 2, (i-1)*N1+ii ) = (ii - (N1+1)/2)*dris;
    end
end

Constant1 = sqrt( dris*dris*lambda*lambda/64/pi^3) ;

H1_los = zeros(Nris,Nt);
for i = 1:Nt
    for ii = 1:Nris
        distance_nt_nm = norm( tx_arr(:,i) - RISPosition(:,ii) );
        
        cos_theta_t_nt_nm = tx_arr(3,i)/distance_nt_nm;
        F = cos_theta_t_nt_nm;
        
        cos_theta_tx_nt_nm = (RISPosition(1,ii) - tx_arr(1,i))/distance_nt_nm;
        F_tx = (cos_theta_tx_nt_nm)^(Gt/2-1);
        
        Amplitude_nt_nm = Constant1*G*Gt*F_tx*F/distance_nt_nm/distance_nt_nm;
        Phase_term_nt_nm = exp( -1j*2*pi*distance_nt_nm/lambda);
        
        H1_los(ii,i) = Amplitude_nt_nm*Phase_term_nt_nm;
    end
end


H2_los = zeros(Nr,Nris);
for i = 1:Nr
    for ii = 1:Nris
        distance_nr_nm = norm( rx_arr(:,i) - RISPosition(:,ii) );
        
        cos_theta_r_nr_nm = rx_arr(3,i)/distance_nr_nm;
        F = cos_theta_r_nr_nm;
        
        cos_theta_rx_nr_nm = (RISPosition(1,ii) - rx_arr(1,i))/distance_nr_nm;
        F_rx = (cos_theta_rx_nr_nm)^(Gr/2-1);
        
        Amplitude_nr_nm = Constant1*G*Gr*F_rx*F/distance_nr_nm/distance_nr_nm;
        Phase_term_nr_nm = exp( -1j*2*pi*distance_nr_nm/lambda);
        
        H2_los(i,ii) = Amplitude_nr_nm*Phase_term_nr_nm;
    end
end





Constant2 = sqrt( lambda*lambda/16/pi^2) ;
Hdir_los = zeros(Nr,Nt);
for i = 1:Nr
    for ii = 1:Nt
        distance_nr_nt = norm( rx_arr(:,i) - tx_arr(:,ii) );
        
        cos_theta_tx_nr_nt = (tx_arr(1,ii) - rx_arr(1,i))/distance_nr_nt;
        F_tx = cos_theta_tx_nr_nt^(Gt/2-1);
        
        cos_theta_rx_nr_nt = (tx_arr(1,ii) - rx_arr(1,i))/distance_nr_nt;
        F_rx = (cos_theta_rx_nr_nt)^(Gr/2-1);
        
        Amplitude_nr_nt = Constant2*Gt*Gr*F_rx*F_tx/distance_nr_nt/distance_nr_nt;
        Phase_term_nr_nt = exp( -1j*2*pi*distance_nr_nt/lambda);
        
        Hdir_los(i,ii) = Amplitude_nr_nt*Phase_term_nr_nt;
    end
end


Hdir = cell(1,no_mat);
H1 = cell(1,no_mat);
H2 = cell(1,no_mat);
for i = 1: no_mat
    Hdir{i} = Hdir_los;
    H1{i} = H1_los;
    H2{i} = H2_los;
end





